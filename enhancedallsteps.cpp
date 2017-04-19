#include "enhancedallsteps.h"

EnhancedAllsteps::EnhancedAllsteps(QObject *parent) : QObject(parent)
{
    secondErrorElevList = new QMap<QString, QPair<double, double> >();
    secondErrorAzimuthList = new QMap<QString, QPair<double, double> >();
    paramRange = new QMap<QString, QPair<double, double> >();
}

EnhancedAllsteps::~EnhancedAllsteps()
{
    delete secondErrorElevList;
    delete secondErrorAzimuthList;
    delete paramRange;
}

void EnhancedAllsteps::slotEnhAllStepsStart()
{
    QString startTime, endTime;
    startTime = QTime::currentTime().toString();
    loadSecondErrorData();

    QVector<QVector2D> hypPoints, sShadows, toPaint;
    QList<QVector3D> p2s;
    QList<double> deg_p2s;
    QVector3D p1, p2, s;

    QImage image;

    double North;
    double R_MIN;
    double deg_p1, deg_p2, err1, err2;

    QString hyp, hypname, picname, ampm;
    int picno;

    int p1_resolution = 30;             // in pixels (30)
    double p2_resolution = 20;          // in deg (20)
    double deltoid_resolution = 9.0;    // (deltoid_points+1)**2 points (9)
    double elev_resolution = 0.2;       // in deg (0.2)
    double azimuth_resolution = 0.1;    // in deg (0.1)
    double second_resolution = 0.25;    // in deg (0.25) .... currently not counted

    QDir folder = QFileDialog::getExistingDirectory();
    QStringList nameList = folder.entryList(QStringList("*_gdeg_masked.tiff"), QDir::Files | QDir::NoDotAndDotDot);
    QStringList imageList;
    QDir::setCurrent(folder.path());
    foreach(QString name, nameList)
        imageList.append(folder.absoluteFilePath(name));

    for (picno = 0; picno < imageList.size(); picno++) { //starting from 0

        picname = imageList.at(picno);

        image = readImage(picname);

        QString imagename = QString(picname.split("/").last()).split(".").first();

        for (int q = 1; q <= 3; q++) { // q = 1..3
            for (int yy = 1; yy <= 4; yy++) { // yy = 1..4

                int stone = q;                      // 1 calcite, 2 cordierite, 3 tourmaline
                if (yy==1) {
                    hyp="/home/denes/Documents/Labor/Viking/ErrorPropagation/l61d0621am.dat";
                    hypname="sol";
                    ampm="am";
                }

                if (yy==2) {
                    hyp="/home/denes/Documents/Labor/Viking/ErrorPropagation/l61d0621pm.dat";
                    hypname="sol";
                    ampm="pm";
                }

                if (yy==3) {
                    hyp="/home/denes/Documents/Labor/Viking/ErrorPropagation/l61d0922am.dat";
                    hypname="equ";
                    ampm="am";
                }

                if (yy==4) {
                    hyp="/home/denes/Documents/Labor/Viking/ErrorPropagation/l61d0922pm.dat";
                    hypname="equ";
                    ampm="pm";
                }

                QVector<QPair<QPair<QVector3D, QVector3D>, QVector3D> > currentFirstErrorPoints, firstErrorPoints_pol;
                QVector<QVector3D> secondErrorPoints_pol, thirdErrorPoints_pol;
                QPair<QVector<QVector2D>, QVector<QVector2D> > sunShadowPoints;
                QVector<int> hist(361);

                usedImage = new QImage(picname);
                usedpixels = 0;

                hypPoints.clear();
                hypPoints = readHyperbolaFromFile(hyp);
                R_MIN = getMinRadius(hypPoints);

                s = detectSun(image);

                double s_elev = (Pi/2 - transform.descartes2Polar(s).y())/Pi*180.0;
                if((yy == 3 || yy == 4) && s_elev > 29)
                    break;

                North = detectNorth(hypPoints, s);
                toPaint.clear();

                int imageWidth = image.width();
                int imageHeight = image.height();

                for(int i=0; i<imageWidth; i+=p1_resolution){
                    for(int j=0; j<imageHeight; j+=p1_resolution){

                        if( QColor(image.pixel(i,j)) != QColor(Qt::red) && QColor(image.pixel(i,j)) != QColor(Qt::green)){
                            usedpixels++;
                            usedImage->setPixelColor(i, j, Qt::blue);
                            QVector2D current_fisheye = transform.draw2Fisheye(QVector2D(i,j),imageWidth);
                            p1 = transform.fisheye2Descartes(current_fisheye);
                            deg_p1 = -(100.0/255.0)*qRed(image.pixel(i,j)) + 100;
                            err1 = getError(deg_p1, stone);
                            QPair<QList<QVector3D>, QList<double> > p2sData = calculatep2s(image, s, p1,45*Pi/180.0, 90*Pi/180.0, p2_resolution*Pi/180.0);
                            p2s = p2sData.first;
                            deg_p2s = p2sData.second;

                            for(int k=0; k<p2s.size(); k++){
                                p2 = p2s.at(k);
                                deg_p2 = deg_p2s.at(k);
                                err2 = getError(deg_p2, stone);
                                QList<QVector3D> sAndps;
                                sAndps.append(s);
                                sAndps.append(p1);
                                sAndps.append(p2);
                                QList<double> firstErrors;
                                firstErrors.append(err1);
                                firstErrors.append(err2);
                                firstErrors.append(deltoid_resolution);
                                currentFirstErrorPoints.clear();
                                currentFirstErrorPoints = getIntersectionsFromFirstStep_pol(sAndps, firstErrors, s);
                                firstErrorPoints_pol.append(currentFirstErrorPoints);
                            }
                        }
                    }
                    QApplication::processEvents();
                    emit signalWriteToList(QString::number((int)((double)i/image.width()*100.0)) + " % of first error calculation ready");
                }
                usedImage->save(picname.remove(".tiff") + "_usedarea.png");
                delete usedImage;
                emit signalWriteToList("First error calculation ready");
                secondErrorPoints_pol = getPointsFromSecondStep_pol(&firstErrorPoints_pol, second_resolution);
                emit signalWriteToList("Second error calculation ready.");
                firstErrorPoints_pol.clear();
                thirdErrorPoints_pol = getPointsFromThirdAndFourthStep_pol(&secondErrorPoints_pol, elev_resolution, azimuth_resolution);
                emit signalWriteToList("Third error calculation ready.");
                secondErrorPoints_pol.clear();
                sunShadowPoints = getSunShadows(&thirdErrorPoints_pol, R_MIN, North);
                sShadows = sunShadowPoints.first;
                toPaint = sunShadowPoints.second;
                emit signalWriteToList("Sun shadow calculation ready.");
                thirdErrorPoints_pol.clear();
                sunShadowPoints.first.clear();
                sunShadowPoints.second.clear();
                hist = calculateNorthErrors(&sShadows, hypPoints, North);
                emit signalWriteToList("North error calculation ready");
                sShadows.clear();
                paint(hypPoints, toPaint, R_MIN, QString(imageList.at(picno)).remove(".tiff") + QString::number(q) + hypname + ampm);
                toPaint.clear();
                writeDatFile(usedpixels, imagename, hypname, ampm, stone, p1_resolution, p2_resolution, deltoid_resolution, second_resolution, elev_resolution, azimuth_resolution, hist, North);

                QApplication::processEvents();
                emit signalWriteToList(QTime::currentTime().toString());
            }
        }

    }

    endTime = QTime::currentTime().toString();
    emit signalWriteToList("Start time: " + startTime);
    emit signalWriteToList("End time: " + endTime);


}

void EnhancedAllsteps::paint(QVector<QVector2D> &hyp_points, QVector<QVector2D> &to_paint, double r_min, QString name)
{
    int size = 1355;
    QPainter painter;
    QPen pen;

    QImage im=QImage(size,size,QImage::Format_ARGB32);
    painter.begin(&im);

    painter.fillRect(0,0,size,size,Qt::white);
    painter.drawEllipse(QPointF(size/2.0,size/2.0),3,3);
    painter.drawEllipse(QPointF(size/2.0,size/2.0),5,5);
    painter.drawEllipse(QPointF(size/2.0,size/2.0),r_min*100, r_min*100);

    pen.setColor(Qt::red);
    painter.setPen(pen);
    int hyp_pointsSize = hyp_points.size();
    for(int q=0; q<hyp_pointsSize; q++)
        painter.drawPoint(size/2.0 + hyp_points.at(q).x()*100, size/2.0 - hyp_points.at(q).y()*100);

    pen.setColor(Qt::black);
    painter.setPen(pen);
    int to_paintSize = to_paint.size();
    for(int q=0; q<to_paintSize; q++)
        painter.drawPoint(size/2.0 + to_paint.at(q).x()*100, size/2.0 - to_paint.at(q).y()*100);

    painter.end();
    im.save(name + ".png");
}

QImage EnhancedAllsteps::readImage(QString filename)
{
    QImage im;
    im = QImage(filename);
    emit signalWriteToList("Image loaded: " + filename);
    return im;
}

double EnhancedAllsteps::getMinRadius(QVector<QVector2D> &hyperbola)
{
    double r_min = 100000.0;
    int hyperbolaSize = hyperbola.size();

    for(int i=0; i< hyperbolaSize; i++){
        qreal hyperbolaLength = hyperbola.at(i).length();
        if (hyperbolaLength < r_min)
            r_min = hyperbolaLength;
    }
    return r_min;
}

QVector<QVector2D> EnhancedAllsteps::readHyperbolaFromFile(QString filename)
{
    QVector<QVector2D> hyp_points;

    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        qDebug("baj");

    QTextStream filestream(&file);
    double x,y;

    while (!filestream.atEnd()){
        QString line = filestream.readLine();
        QTextStream stream(&line);
        stream >> x >> y;
        hyp_points.append(QVector2D(x,y));
    }

    file.close();
    emit signalWriteToList(filename + QString(" read ready"));
    return hyp_points;
}

double EnhancedAllsteps::getError(double deg_of_pol, int stein)
{
    if (stein == 1) {                                             // calcite 2 contrast
        return 1 / ( (0.00652452) * deg_of_pol + 0.0175308 ) + 4.71625;
    }

    if (stein == 2) {                                              // cordierite
        return 1 / ( (0.00422263) * deg_of_pol + 0.0117355 ) + 0.55487;
    }

    if (stein == 3)                                                // tourmaline
        return 1 / ( (0.00266655) * deg_of_pol + 0.0102908 ) - 0.669271;

    else{
        emit signalWriteToList("Not valid stone!");
        return 0;
    }
}

void EnhancedAllsteps::writeDatFile(int usedpixels, QString imname, QString hyperbola, QString timeOfDay, int stein, int p1_res, double p2_res, double deltoid_res, double second_error_res, double elev_res, double azimuth_res, QVector<int> &histogram, double north)
{
    QString s;
    if(stein == 1)
        s = QString("cal");
    if(stein == 2)
        s = QString("cord");
    if(stein == 3)
        s = QString("tour");

    QString filename = imname + "_" + hyperbola + "_" + timeOfDay + "_" + s + "_" + QString::number(p1_res) + "px_" + QString::number(p2_res) + "deg_" + QString::number((deltoid_res+1)*(deltoid_res+1)) + "pts_" /*+ QString::number(second_error_res) + "deg_"*/ + QString::number(elev_res) + "deg_" + QString::number(azimuth_res) + "deg";
    QFile file(filename + ".csv");
    if (!file.open(QIODevice::WriteOnly))
        qDebug("baj");
    QTextStream stream(&file);

    stream << "#N = " << histogram.at(360) << "\n";
    stream << "#North = " << north/Pi*180.0 << " deg\n";
    stream << "#Used px ratio = " << QString::number(usedpixels/980100.0) << "\n"; //990*990

    int histSize = histogram.size();
    for (int i=0; i<histSize-1; i++)
        stream << i-180 << "\t" << histogram.at(i) << "\n";

    file.close();
    emit signalWriteToList(filename + QString(" ready"));
}

QVector3D EnhancedAllsteps::detectSun(QImage &im)
{
    QVector3D sun;
    int imageWidth = im.width();
    int imageHeight = im.height();

    for(int i=0; i<imageWidth; i++){
        for (int j=0; j<imageHeight; j++){
            if (QColor(im.pixel(i,j)) == QColor(Qt::green)){
                QVector2D current_fisheye = transform.draw2Fisheye(QVector2D(i,j),imageWidth);
                sun = transform.fisheye2Descartes(current_fisheye);
                emit signalWriteToList("sun detected " + QString::number(i) + " " + QString::number(j));
                emit signalWriteToList(QString::number((Pi/2 - transform.descartes2Polar(sun).y())/Pi*180.0));
            }
        }
    }
    return sun;
}

double EnhancedAllsteps::detectNorth(QVector<QVector2D> &hyp_points, QVector3D &sun)
{
    double north;
    QVector<double> nerrors;
    QVector<QVector2D> minimums;
    QVector2D INTERSECT;

    QVector3D sun_pol = transform.descartes2Polar(sun);
    QVector2D sun_shadow = QVector2D(tan(sun_pol.y())*cos(-sun_pol.z()) , tan(sun_pol.y())*sin(-sun_pol.z()));

    minimums.clear();
    nerrors.clear();

    int hyp_pointsSize = hyp_points.size();

    for (int j=1; j<hyp_pointsSize-1; j++){
        if( (fabs(sun_shadow.length()-hyp_points.at(j).length())<fabs(sun_shadow.length()-hyp_points.at(j-1).length()))&&(fabs(sun_shadow.length()-hyp_points.at(j).length())<fabs(sun_shadow.length()-hyp_points.at(j+1).length()))){
            minimums.append(hyp_points.at(j));
            nerrors.append(acos(QVector2D::dotProduct(hyp_points.at(j).normalized(), sun_shadow.normalized())));
        }
    }

    INTERSECT = minimums.at(0);
    north = nerrors.at(0);

    int minimumsSize = minimums.size();
    for(int j=0; j<minimumsSize; j++)
        if(nerrors.at(j) < north){
            INTERSECT = minimums.at(j);
            north = nerrors.at(j);
        }


    if (QVector3D::crossProduct(QVector3D (sun_shadow.x(), sun_shadow.y(),0.0), QVector3D (INTERSECT.x(), INTERSECT.y(), 0.0) ).z() < 0)
        north *= (-1.0);

    emit signalWriteToList("north detected: " + QString::number(north*180/Pi) + " deg");
    return north;
}

QPair<QList<QVector3D>, QList<double> > EnhancedAllsteps::calculatep2s(QImage &im, QVector3D &sun, QVector3D &p1, double min, double max, double res)
{
    QPair<QList<QVector3D>, QList<double> > p2s_data;
    QVector3D axis,v;
    axis = QVector3D::crossProduct(p1,sun).normalized();
    int k,l;

    int imageWidth = im.width();

    for(double i=min; i<=max; i += res){
        for(double j=0.0; j<2*Pi; j += res){
            v = transform.rotate ( transform.rotate(p1,axis,i) , p1.normalized() , j);
            if(v.z()>0){
                k = (int)transform.fisheye2Draw(transform.descartes2Fisheye(v),imageWidth).x();
                l = (int)transform.fisheye2Draw(transform.descartes2Fisheye(v),imageWidth).y();
                if(k>=0 && k<imageWidth && l>=0 && l<imageWidth && QColor(im.pixel(k,l))!=QColor(Qt::red)){
                    usedpixels++;
                    usedImage->setPixelColor(k, l, Qt::blue);
                    p2s_data.first.append(v);
                    p2s_data.second.append(-(100.0/255.0)*qRed(im.pixel(k,l)) + 100);
                }
            }
        }
    }
    return p2s_data;
}

QVector3D EnhancedAllsteps::intersectionOfGreatCircles(QVector3D &GC1A, QVector3D &GC1B, QVector3D &GC2A, QVector3D &GC2B, QVector3D &sun)
{
    QVector3D D,result;
    QVector3D axis1 = QVector3D::crossProduct(GC1A, GC1B).normalized();
    QVector3D axis2 = QVector3D::crossProduct(GC2A, GC2B).normalized();

    if (axis1==axis2){
        result = sun;
    }
    else{
        D = (QVector3D::crossProduct(axis1, axis2)).normalized();
        if (D.z() >= 0)
            result = D;
        else
            result = -D;
    }
    return result;
}

QVector<int> EnhancedAllsteps::calculateNorthErrors(QVector<QVector2D> *sun_shadows, QVector<QVector2D> &hyp_points, double north)
{
    QVector<int> histogram(361);
    QVector<double> nerrors;
    double NERROR;
    QVector<QVector2D> minimums;
    QVector2D INTERSECT;
    int h;

    int sShadowsSize = sun_shadows->size();
    int hyp_pointsSize = hyp_points.size();

    for(int i=0; i<sShadowsSize; i++){
        minimums.clear();
        nerrors.clear();

        QVector2D current_shadow = sun_shadows->at(i);
        sun_shadows->replace(i, transform.rotate2D(current_shadow, north));

        for (int j=1; j<hyp_pointsSize-1; j++){
            if( (fabs(sun_shadows->at(i).length()-hyp_points.at(j).length())<fabs(sun_shadows->at(i).length()-hyp_points.at(j-1).length()))&&(fabs(sun_shadows->at(i).length()-hyp_points.at(j).length())<fabs(sun_shadows->at(i).length()-hyp_points.at(j+1).length()))){
                minimums.append(hyp_points.at(j));
                nerrors.append(acos(QVector2D::dotProduct(hyp_points.at(j).normalized(), sun_shadows->at(i).normalized())));
            }
        }

        if (minimums.size() > 0){
            INTERSECT = minimums.at(0);
            NERROR = nerrors.at(0);

            int minimumsSize = minimums.size();
            for(int j=0; j<minimumsSize; j++)
                if(nerrors.at(j) < NERROR){
                    INTERSECT = minimums.at(j);
                    NERROR = nerrors.at(j);
                }

            if (QVector3D::crossProduct(sun_shadows->at(i).toVector3D(), INTERSECT.toVector3D()).z() < 0)
                NERROR *= (-1.0);

            NERROR = NERROR/Pi*180.0;
            h = (int)(NERROR + 180) % 360;
            if(h >= 0){
                histogram[h]++;
                histogram[360]++;
            }
        }
        QApplication::processEvents();
        if(i % 100000 == 0)
            emit signalWriteToList("North error calculation " + QString::number(100*(double)i/(double)sShadowsSize) + " % ready.");
    }
    return histogram;
}

QVector<QPair<QPair<QVector3D, QVector3D>, QVector3D> > EnhancedAllsteps::getIntersectionsFromFirstStep_pol(QList<QVector3D> &sunAndPoints, QList<double> &firstStepErrors, QVector3D &sun)
{
    QVector<QPair<QPair<QVector3D, QVector3D>, QVector3D> > v_pol; //first: point1, point2 second: v_pol
    QPair<QPair<QVector3D, QVector3D>, QVector3D> current_v_pol;
    current_v_pol.first.first = sunAndPoints.at(1);
    current_v_pol.first.second = sunAndPoints.at(2);
    if(sunAndPoints.size() == 3 && firstStepErrors.size() == 3){
        QVector3D GC1_point, GC2_point, v;
        QVector3D p1, p2;
        QVector3D p1Plus, p1Minus, n_p1Plus, n_p1Minus;
        QVector3D p2Plus, p2Minus, n_p2Plus, n_p2Minus;

        p1Plus  = transform.rotate (sunAndPoints.at(1), transform.rotate(QVector3D::crossProduct(sunAndPoints.at(1), sunAndPoints.at(0)).normalized(), sunAndPoints.at(1),  firstStepErrors.at(0)*Pi/180.0), Pi/10.0 );
        p1Minus = transform.rotate (sunAndPoints.at(1), transform.rotate(QVector3D::crossProduct(sunAndPoints.at(1), sunAndPoints.at(0)).normalized(), sunAndPoints.at(1), -firstStepErrors.at(0)*Pi/180.0), Pi/10.0 );
        p2Plus  = transform.rotate (sunAndPoints.at(2), transform.rotate(QVector3D::crossProduct(sunAndPoints.at(2), sunAndPoints.at(0)).normalized(), sunAndPoints.at(2),  firstStepErrors.at(1)*Pi/180.0), Pi/10.0 );
        p2Minus = transform.rotate (sunAndPoints.at(2), transform.rotate(QVector3D::crossProduct(sunAndPoints.at(2), sunAndPoints.at(0)).normalized(), sunAndPoints.at(2), -firstStepErrors.at(1)*Pi/180.0), Pi/10.0 );

        n_p1Plus  = QVector3D::crossProduct(sunAndPoints.at(1), p1Plus).normalized();
        n_p1Minus = QVector3D::crossProduct(sunAndPoints.at(1), p1Minus).normalized();
        n_p2Plus  = QVector3D::crossProduct(sunAndPoints.at(2), p2Plus).normalized();
        n_p2Minus = QVector3D::crossProduct(sunAndPoints.at(2), p2Minus).normalized();

        for(double i=0.0; i<=firstStepErrors.at(2); i+=1){
            GC1_point = transform.rotate (sunAndPoints.at(1), transform.rotate(QVector3D::crossProduct(sunAndPoints.at(1), sunAndPoints.at(0)).normalized(), sunAndPoints.at(1), (-firstStepErrors.at(0)+2.0*i*firstStepErrors.at(0)/firstStepErrors.at(2))*Pi/180.0), Pi/10.0 );

            for (double j=0.0; j<=firstStepErrors.at(2); j+=1){

                GC2_point = transform.rotate (sunAndPoints.at(2), transform.rotate(QVector3D::crossProduct(sunAndPoints.at(2), sunAndPoints.at(0)).normalized(), sunAndPoints.at(2), (-firstStepErrors.at(1)+2.0*j*firstStepErrors.at(1)/firstStepErrors.at(2))*Pi/180.0), Pi/10.0 );
                p1 = sunAndPoints.at(1);
                p2 = sunAndPoints.at(2);
                v = intersectionOfGreatCircles(p1, GC1_point, p2, GC2_point, sun);
                if(QVector3D::dotProduct(n_p1Plus,sunAndPoints.at(0))*QVector3D::dotProduct(n_p1Plus,v) > 0 && QVector3D::dotProduct(n_p1Minus,sunAndPoints.at(0))*QVector3D::dotProduct(n_p1Minus,v) > 0 && QVector3D::dotProduct(n_p2Plus,sunAndPoints.at(0))*QVector3D::dotProduct(n_p2Plus,v) > 0 && QVector3D::dotProduct(n_p2Minus,sunAndPoints.at(0))*QVector3D::dotProduct(n_p2Minus,v) > 0){
                    current_v_pol.second = transform.descartes2Polar(v);
                    v_pol.append(current_v_pol);
                }
            }
        }

    }
    else{
        emit signalWriteToList("Error in getIntersectionsFromFirstStep_pol: condition not fulfilled.");
        current_v_pol.second = QVector3D(-999,-999,-999);
        v_pol.append(current_v_pol);
        return v_pol;
    }
    return v_pol;
}

QVector<QVector3D> EnhancedAllsteps::getPointsFromSecondStep_pol(QVector<QPair<QPair<QVector3D, QVector3D>, QVector3D> > *intersectPoints_pol, double seconderror_resolution)
{
    QVector<QVector3D> secondStepPoints;
    if(!intersectPoints_pol->isEmpty()){
        int num = intersectPoints_pol->size();
        QVector3D v_des;
        QPair<QPair<QVector3D, QVector3D>, QVector3D> currentPoint;
        QList<double> parameters;
        QVector3D new_v_pol;
        QVector3D firstSecond_sun;
        QPair<double, double> secondElevErrors, secondAzimuthErrors;
        for(int i = 0; i < num; i++){
            currentPoint = intersectPoints_pol->at(i);
            double v_elev_deg = (Pi/2 - currentPoint.second.y())/Pi*180.0;
            v_des = transform.polar2Descartes(currentPoint.second);
            double delta_deg = acos(QVector3D::dotProduct(QVector3D::crossProduct(currentPoint.first.first, v_des).normalized(), QVector3D::crossProduct(currentPoint.first.second, v_des).normalized())) * 180.0/Pi;
            double th1_deg = acos(QVector3D::dotProduct(currentPoint.first.first.normalized(), v_des.normalized())) * 180.0/Pi;
            double th2_deg = acos(QVector3D::dotProduct(currentPoint.first.second.normalized(), v_des.normalized())) * 180.0/Pi;

            parameters.clear();
            parameters.append(delta_deg);
            parameters.append(v_elev_deg);
            parameters.append(th1_deg);
            parameters.append(th2_deg);

            secondElevErrors = getSeconStepError(parameters, secondErrorElevList, paramRange, secondErrorElevListKeysNum);
            secondAzimuthErrors = getSeconStepError(parameters, secondErrorAzimuthList, paramRange, secondErrorAzimuthListKeysNum);

            if((int)secondElevErrors.first != -999 && (int)secondAzimuthErrors.first != -999){

                new_v_pol.setX(currentPoint.second.x());
                new_v_pol.setY(currentPoint.second.y() - (secondElevErrors.first * Pi/180.0)); /*theta measured from vertical*/
                new_v_pol.setZ(currentPoint.second.z() + (secondAzimuthErrors.first * Pi/180.0));
                firstSecond_sun.setX(new_v_pol.x());
/*
                for(double el = -secondElevErrors.second; el <= secondElevErrors.second; el += seconderror_resolution){
                    for(double az = -secondAzimuthErrors.second; az < secondAzimuthErrors.second; az += seconderror_resolution){
                        firstSecond_sun.setY(new_v_pol.y() + (el * Pi/180.0));
                        firstSecond_sun.setZ(new_v_pol.z() + (az * Pi/180.0));
                        secondStepPoints.append(firstSecond_sun);
                    }
                }
*/
                firstSecond_sun.setY(new_v_pol.y());
                firstSecond_sun.setZ(new_v_pol.z());
                secondStepPoints.append(firstSecond_sun);
            }
            QApplication::processEvents();
            if(i % 100000 == 0)
                emit signalWriteToList("Second step " + QString::number(100*(double)i/(double)num) + " % ready.");
        }
    }
    else{
        emit signalWriteToList("Error in getPointsFromSecondStep: condition not fulfilled.");
        secondStepPoints.append(QVector3D(-999,-999,-999));
        return secondStepPoints;
    }
    return secondStepPoints;
}

QVector<QVector3D> EnhancedAllsteps::getPointsFromThirdAndFourthStep_pol(QVector<QVector3D> *secondErrorPoint_pol, double elev_res_deg, double azimuth_res_deg)
{
    QVector<QVector3D> thirdStepPoints;
    if(!secondErrorPoint_pol->isEmpty()){
        double count=0;
        foreach(QVector3D currentpoint, *secondErrorPoint_pol){
            double thirdError;
            double firstSecond_elev_deg = (Pi/2.0 - currentpoint.y()) * 180.0/Pi;
            thirdError = 0.112482*pow(firstSecond_elev_deg, 0.784039); //best
            if(thirdError < 0)
                thirdError = 0;

            QVector3D sunWithErrors;
            sunWithErrors.setX(currentpoint.x());

            for(double elev = -thirdError; elev <= thirdError; elev += elev_res_deg){
                sunWithErrors.setY(currentpoint.y() + elev*Pi/180.0);
                if((Pi/2.0 - sunWithErrors.y()) < 0) {
                    elev = -currentpoint.y()*180.0/Pi;
                    sunWithErrors.setY(currentpoint.y() + elev*Pi/180.0);
                }
                double azimuthError_deg = getAzimuthError(sunWithErrors.y()*180.0/Pi);
                for(double azError = -azimuthError_deg; azError <= azimuthError_deg; azError += azimuth_res_deg){
                    sunWithErrors.setZ(currentpoint.z() + azError*Pi/180.0);
                    thirdStepPoints.append(sunWithErrors);
                }
            }
            count++;
            QApplication::processEvents();
            if((int)count % 100000 == 0)
                emit signalWriteToList("Third step " + QString::number(100*count/(double)secondErrorPoint_pol->size()) + " % ready.");
        }
    }
    else{
        emit signalWriteToList("Error in getPointsFromThirdAndFourthStep_pol: condition not fulfilled.");
        thirdStepPoints.append(QVector3D(-999,-999,-999));
        return thirdStepPoints;
    }
    return thirdStepPoints;
}

QPair<QVector<QVector2D>, QVector<QVector2D> > EnhancedAllsteps::getSunShadows(QVector<QVector3D> *thirdErrorPoints_pol, double r_min, double north)
{
    QPair<QVector<QVector2D>, QVector<QVector2D> > sunShadows;
    double count = 0;
    if(!thirdErrorPoints_pol->isEmpty()){
        QVector2D projected, paint;
        foreach(QVector3D currentpoint, *thirdErrorPoints_pol){
            projected = QVector2D( tan(currentpoint.y())*cos(-currentpoint.z()) , tan(currentpoint.y())*sin(-currentpoint.z()) );
            if(currentpoint.y()<Pi/2.0-Pi/16.0  && projected.length()>r_min){ // horizonthoz ne legyen tul kozel (nagyon hosszu arnyek)
                sunShadows.first.append(projected);
                paint = sunShadows.first.last();
                paint = transform.rotate2D(paint, north);
                sunShadows.second.append(paint);
            }
            count++;
            QApplication::processEvents();
            if((int)count % 100000 == 0)
                emit signalWriteToList("Sunshadow " + QString::number(100*count/(double)thirdErrorPoints_pol->size()) + " % ready.");
        }
    }
    else{
        emit signalWriteToList("Error in getSunShadows: condition not fulfilled.");
        sunShadows.first.append(QVector2D(-999,-999));
        sunShadows.second.append(QVector2D(-999,-999));
        return sunShadows;
    }
    return sunShadows;
}

double EnhancedAllsteps::getAzimuthError(double elevation)
{
    double result = 0.000238133*elevation*elevation + 1.1829;
}

void EnhancedAllsteps::loadSecondErrorData()
{
    secondErrorElevList->clear();
    secondErrorAzimuthList->clear();
    paramRange->clear();

    QFile file("/home/denes/Documents/Labor/Viking/ErrorPropagation/result_noantisun_ranges.csv");
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        qDebug("cannot open");
    QTextStream stream(&file);
    QString str;
    while(!stream.atEnd()){
        QString line = stream.readLine();
        if(!line.startsWith("#")){
            QStringList lineList = line.split("\t");
            str = lineList.at(0);
            secondErrorElevList->insert(str, QPair<double, double> (QString(lineList.at(1)).toDouble(), QString(lineList.at(2)).toDouble()));
            secondErrorAzimuthList->insert(str, QPair<double, double> (QString(lineList.at(3)).toDouble(), QString(lineList.at(4)).toDouble()));
        }
    }

    for(int j = 0; j < secondErrorAzimuthList->keys().size(); j++){
        QString currentkey = secondErrorAzimuthList->keys().at(j);
        QStringList keyparts = currentkey.split(" ");

        for(int i = 0; i < keyparts.size(); i++) {
            QStringList currentPart = QString(keyparts.at(i)).split("-");
            double min = QString(currentPart.takeFirst()).toDouble();
            double max = QString(currentPart.takeFirst()).toDouble();
            if(i == 0){
                paramRange->insert("delta_" + QString::number(j), QPair<double, double> (min, max));
            }
            if(i == 1){
                paramRange->insert("elev_" + QString::number(j), QPair<double, double> (min, max));
            }
            if(i == 2){
                paramRange->insert("th1_" + QString::number(j), QPair<double, double> (min, max));
            }
            if(i == 3){
                paramRange->insert("th2_" + QString::number(j), QPair<double, double> (min, max));
            }
        }
    }

    secondErrorElevListKeysNum = secondErrorElevList->keys().size();
    secondErrorAzimuthListKeysNum = secondErrorAzimuthList->keys().size();

    emit signalWriteToList("Second error data loaded.");
}

QPair<double, double> EnhancedAllsteps::getSeconStepError(QList<double> &params, QMap<QString, QPair<double, double> > *secondErrorMap, QMap<QString, QPair<double, double> > *rangeMap, double keysNum)
{
    QPair<double, double> seconderror(-999, -999); /*to be easy to recognise*/
    if(!params.isEmpty() && params.size() == 4){
        QString currentkey;
        for(int j = 0; j < keysNum; j++){
            currentkey = secondErrorMap->keys().at(j);
            if(params.at(0) >= rangeMap->value("delta_" + QString::number(j)).first && params.at(0) <= rangeMap->value("delta_" + QString::number(j)).second && params.at(1) >= rangeMap->value("elev_" + QString::number(j)).first && params.at(1) <= rangeMap->value("elev_" + QString::number(j)).second && params.at(2) >= rangeMap->value("th1_" + QString::number(j)).first && params.at(2) <= rangeMap->value("th1_" + QString::number(j)).second && params.at(3) >= rangeMap->value("th2_" + QString::number(j)).first && params.at(3) <= rangeMap->value("th2_" + QString::number(j)).second){
                seconderror.first = secondErrorMap->value(currentkey).first;
                seconderror.second = secondErrorMap->value(currentkey).second;
                break;
            }
        }
    }
    return seconderror;
}
