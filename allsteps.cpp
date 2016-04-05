#include "allsteps.h"

AllSteps::AllSteps(QObject *parent) : QObject(parent)
{
    image = new QImage();
}

void AllSteps::slotAllStepsStart()
{
    emit signalWriteToList(QTime::currentTime().toString());
    loadSecondErrorData();

    QString hyp, hypname, picname, ampm, bw;
    int picno;

    int p1_resolution = 20;             // in pixels (20)
    double p2_resolution = 15;          // in deg (15)
    double deltoid_resolution = 19.0;    // (deltoid_points+1)**2 points (19)
    double elev_resolution = 0.2;       // in deg (0.2)

    QDir folder = QFileDialog::getExistingDirectory();
    QStringList nameList = folder.entryList(QStringList(), QDir::Files | QDir::NoDotAndDotDot);
    QStringList imageList;
    foreach(QString name, nameList)
        imageList.append(folder.absoluteFilePath(name));


    for (picno=0; picno < 1/*imageList.size()*/; picno++) {

        picname = imageList.at(picno);

        readImage(picname);
        QString imagename = QString(picname.split("/").last()).split(".").first();

        for (int m = 1; m <= 1; m++){ // m = 1..2
            for (int q = 1; q <= 1; q++) { // q = 1..4
                for (int yy = 1; yy <= 1; yy++) { // yy = 1..4
                    for (int z = 1; z <= 1; z++){ // z = 1..3
                        if (z == 1) bw = "best";
                        if (z == 2) bw = "second";
                        if (z == 3) bw = "worst";

                        int stone = q;                      // 1 B, 2 C, 3 D, 4 E
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

                        readHyperbolaFromFile(QString(hyp));
                        R_MIN = getMinRadius(hypPoints);

                        detectSun();
                        detectNorth();
                        initializeHist();
                        toPaint.clear();

                        int imageWidth = image->width();
                        int imageHeight = image->height();

                        for(int i=0; i<imageWidth; i+=p1_resolution){
                            for(int j=0; j<imageHeight; j+=p1_resolution){

                                if( QColor(image->pixel(i,j)) != QColor(Qt::red) && QColor(image->pixel(i,j)) != QColor(Qt::green)){
                                    QVector2D current_fisheye = transform.draw2Fisheye(QVector2D(i,j),imageWidth);
                                    p1 = transform.fisheye2Descartes(current_fisheye);

                                    deg_p1 = -(100.0/255.0)*qRed(image->pixel(i,j)) + 100;
                                    err1 = getError(m, deg_p1, stone);
                                    calculatep2s(p1,45*Pi/180.0, 90*Pi/180.0, p2_resolution*Pi/180.0);

                                    for(int k=0; k<p2s.size(); k++){
                                        p2 = p2s.at(k);
                                        deg_p2 = deg_p2s.at(k);
                                        err2 = getError(m, deg_p2, stone);
                                        calculateSunPositions(s,p1,p2,err1,err2,deltoid_resolution, elev_resolution, bw);
                                        calculateNorthErrors();
                                        QCoreApplication::processEvents();

                                    }
                                }
                            }
                            QApplication::processEvents();
                            emit signalWriteToList(QString::number((int)((double)i/image->width()*100.0)) + " % ready");
                        }
                        paint();
                        QApplication::processEvents();
                        emit signalWriteToList("100 % ready\n");
                        emit signalWriteToList(QTime::currentTime().toString());

                        writeDatFile(imagename, hypname, ampm, stone, p1_resolution, p2_resolution, deltoid_resolution, m, bw);
                    }
                }

            }
        }
        delete image;
    }
}

double AllSteps::getMinRadius(QVector<QVector2D> hyperbola)
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

void AllSteps::initializeHist()
{
    hist.resize(361);
    int histSize = hist.size();
    for (int i=0; i<histSize; i++)
        hist[i] = 0;
}

double AllSteps::getError(int mode, double deg_of_pol, int stein)
{
    if (mode == 1){                                             //equal
        if (stein == 1) {                                             // calcite B
            if (deg_of_pol < 36.9)
                return 22.5;
            else
                return 1 / ( (0.000421467) * deg_of_pol + 0.0131553 ) - 11.3618;
        }

        if (stein == 2) {                                              // calcite C
            if (deg_of_pol < 7.1)
                return 22.5;
            else
                return 1 / ( (0.0127314) * deg_of_pol - 0.0626234 ) + 6.49337;
        }

        if (stein == 3)                                                // calcite D
            return 1 / ( (0.0112355) * deg_of_pol - 0.3914 ) + 6.53672;

        if (stein == 4)                                                // calcite E
            return 1 / ( (0.00294135) * deg_of_pol + 0.0253163 ) + 2.57876;

        else{
            emit signalWriteToList("Not valid stone!");
            return 0;
        }
    }

    if (mode == 2){                                                 //contrast
        if (stein == 1) {                                             // calcite B

            return 1 / ( (0.00106223) * deg_of_pol + 0.0283133 ) - 3.04705;
        }

        if (stein == 2) {                                              // calcite C

            return 1 / ( (0.00652452) * deg_of_pol + 0.0175308 ) + 4.71625;
        }

        if (stein == 3)                                                // calcite D
            return 1 / ( (0.0292764) * deg_of_pol - 0.0585575 ) + 13.301;

        if (stein == 4)                                                // calcite E
            return 1 / ( (0.00251985) * deg_of_pol + 0.0350572 ) + 1.63867;

        else{
            emit signalWriteToList("Not valid stone!");
            return 0;
        }
    }
    else
        return 0;
}

void AllSteps::readHyperbolaFromFile(QString filename)
{
    hypPoints.clear();

    QFile file(filename);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        qDebug("baj");

    QTextStream filestream(&file);
    double x,y;

    while (!filestream.atEnd()){
        QString line = filestream.readLine();
        QTextStream stream(&line);
        stream >> x >> y;
        hypPoints.append(QVector2D(x,y));
    }

    file.close();
    emit signalWriteToList(filename + QString(" read ready"));
}

void AllSteps::writeDatFile(QString imname, QString hyperbola, QString timeOfDay, int stein, int p1_res, double p2_res, double deltoid_res, int mode, QString bestworst)
{
    QString s, Mode;
    if(stein == 1)
        s = QString("B");
    if(stein == 2)
        s = QString("C");
    if(stein == 3)
        s = QString("D");
    if(stein == 4)
        s = QString("E");

    if (mode == 1)
        Mode = "equal";
    if (mode == 2)
        Mode = "contrast";


    QString filename = imname + "_" + hyperbola + "_" + timeOfDay + "_" + s + "_" + QString::number(p1_res) + "px_" + QString::number(p2_res) + "deg_" + QString::number((deltoid_res+1)*(deltoid_res+1)) + "pts_" + Mode + "_" + bestworst;
    QFile file(filename + ".csv");
    if (!file.open(QIODevice::WriteOnly))
        qDebug("baj");
    QTextStream stream(&file);

    stream << "#N = " << hist.at(360) << "\n";
    stream << "#North = " << North/Pi*180.0 << " deg\n";

    int histSize = hist.size();
    for (int i=0; i<histSize-1; i++)
        stream << i-180 << "\t" << hist.at(i) << "\n";

    file.close();
    emit signalWriteToList(filename + QString(" ready"));
}

void AllSteps::readImage(QString filename)
{
    image = new QImage(filename);
    emit signalWriteToList("Image loaded: " + filename);
}

void AllSteps::detectSun()
{
    int imageWidth = image->width();
    int imageHeight = image->height();

    for(int i=0; i<imageWidth; i++){
        for (int j=0; j<imageHeight; j++){
            if (QColor(image->pixel(i,j)) == QColor(Qt::green)){
                QVector2D current_fisheye = transform.draw2Fisheye(QVector2D(i,j),imageWidth);
                s = transform.fisheye2Descartes(current_fisheye);
                emit signalWriteToList("sun detected " + QString::number(i) + " " + QString::number(j));
                emit signalWriteToList(QString::number((Pi/2 - transform.descartes2Polar(s).y())/Pi*180.0));
            }
        }
    }
}

void AllSteps::detectNorth()
{
    QVector<double> nerrors;
    QVector<QVector2D> minimums;
    QVector2D INTERSECT;

    QVector3D sun_pol = transform.descartes2Polar(s);
    QVector2D sun_shadow = QVector2D(tan(sun_pol.y())*cos(-sun_pol.z()) , tan(sun_pol.y())*sin(-sun_pol.z()));

    minimums.clear();
    nerrors.clear();

    int hypPointsSize = hypPoints.size();

    for (int j=1; j<hypPointsSize-1; j++){
        if( (fabs(sun_shadow.length()-hypPoints.at(j).length())<fabs(sun_shadow.length()-hypPoints.at(j-1).length()))&&(fabs(sun_shadow.length()-hypPoints.at(j).length())<fabs(sun_shadow.length()-hypPoints.at(j+1).length()))){
            minimums.append(hypPoints.at(j));
            nerrors.append(acos(QVector2D::dotProduct(hypPoints.at(j).normalized(), sun_shadow.normalized())));
        }
    }

    INTERSECT = minimums.at(0);
    North = nerrors.at(0);

    int minimumsSize = minimums.size();
    for(int j=0; j<minimumsSize; j++)
        if(nerrors.at(j) < North){
            INTERSECT = minimums.at(j);
            North = nerrors.at(j);
        }


    if (QVector3D::crossProduct(QVector3D (sun_shadow.x(), sun_shadow.y(),0.0), QVector3D (INTERSECT.x(), INTERSECT.y(), 0.0) ).z() < 0)
        North *= (-1.0);

    emit signalWriteToList("north detected: " + QString::number(North*180/Pi) + " deg");
}

void AllSteps::calculatep2s(QVector3D p1, double min, double max, double res)
{
    p2s.clear();
    deg_p2s.clear();
    QVector3D axis,v;
    axis = QVector3D::crossProduct(p1,s).normalized();
    int k,l;

    int imageWidth = image->width();

    for(double i=min; i<=max; i += res){
        for(double j=0.0; j<2*Pi; j += res){
            v = transform.rotate ( transform.rotate(p1,axis,i) , p1.normalized() , j);
            if(v.z()>0){
                k = (int)transform.fisheye2Draw(transform.descartes2Fisheye(v),imageWidth).x();
                l = (int)transform.fisheye2Draw(transform.descartes2Fisheye(v),imageWidth).y();
                if(k>=0 && k<imageWidth && l>=0 && l<imageWidth && QColor(image->pixel(k,l))!=QColor(Qt::red)){
                     p2s.append(v);
                     deg_p2s.append(-(100.0/255.0)*qRed(image->pixel(k,l)) + 100);
                }
            }
        }
    }
}

QVector3D AllSteps::intersectionOfGreatCircles(QVector3D GC1A, QVector3D GC1B, QVector3D GC2A, QVector3D GC2B)
{
    QVector3D D,result;
    QVector3D axis1 = QVector3D::crossProduct(GC1A, GC1B).normalized();
    QVector3D axis2 = QVector3D::crossProduct(GC2A, GC2B).normalized();

    if (axis1==axis2){
        result = s;
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

void AllSteps::calculateSunPositions(QVector3D sun, QVector3D point1, QVector3D point2, double error1, double error2, double res, double elev_res_deg, QString bw)
{
    QVector2D projected;
    QVector2D paint;

    sShadows.clear();
    QVector3D GC1_point, GC2_point, v;
    QVector3D p1Plus, p1Minus, n_p1Plus, n_p1Minus;
    QVector3D p2Plus, p2Minus, n_p2Plus, n_p2Minus;

    p1Plus  = transform.rotate (point1, transform.rotate(QVector3D::crossProduct(point1, sun).normalized(), point1,  error1*Pi/180.0), Pi/10.0 );
    p1Minus = transform.rotate (point1, transform.rotate(QVector3D::crossProduct(point1, sun).normalized(), point1, -error1*Pi/180.0), Pi/10.0 );
    p2Plus  = transform.rotate (point2, transform.rotate(QVector3D::crossProduct(point2, sun).normalized(), point2,  error2*Pi/180.0), Pi/10.0 );
    p2Minus = transform.rotate (point2, transform.rotate(QVector3D::crossProduct(point2, sun).normalized(), point2, -error2*Pi/180.0), Pi/10.0 );

    n_p1Plus  = QVector3D::crossProduct(point1, p1Plus).normalized();
    n_p1Minus = QVector3D::crossProduct(point1, p1Minus).normalized();
    n_p2Plus  = QVector3D::crossProduct(point2, p2Plus).normalized();
    n_p2Minus = QVector3D::crossProduct(point2, p2Minus).normalized();

    for(double i=0.0; i<=res; i+=1){
        for (double j=0.0; j<=res; j+=1){

            GC1_point = transform.rotate (point1, transform.rotate(QVector3D::crossProduct(point1, sun).normalized(), point1, (-error1+2.0*i*error1/res)*Pi/180.0), Pi/10.0 );
            GC2_point = transform.rotate (point2, transform.rotate(QVector3D::crossProduct(point2, sun).normalized(), point2, (-error2+2.0*j*error2/res)*Pi/180.0), Pi/10.0 );

            v = intersectionOfGreatCircles(point1, GC1_point, point2, GC2_point);

            if(QVector3D::dotProduct(n_p1Plus,sun)*QVector3D::dotProduct(n_p1Plus,v) > 0 && QVector3D::dotProduct(n_p1Minus,sun)*QVector3D::dotProduct(n_p1Minus,v) > 0 && QVector3D::dotProduct(n_p2Plus,sun)*QVector3D::dotProduct(n_p2Plus,v) > 0 && QVector3D::dotProduct(n_p2Minus,sun)*QVector3D::dotProduct(n_p2Minus,v) > 0){

                v = transform.descartes2Polar(v);

                /*second step error*/
                double v_elev_deg = (Pi/2 - v.y())/Pi*180.0;
                QVector3D v_des = transform.polar2Descartes(v);
                double delta_deg = acos(QVector3D::dotProduct(QVector3D::crossProduct(point1, v_des).normalized(), QVector3D::crossProduct(point2, v_des).normalized())) * 180.0/Pi;
                double th1_deg = acos(QVector3D::dotProduct(point1.normalized(), v_des.normalized())) * 180.0/Pi;
                double th2_deg = acos(QVector3D::dotProduct(point2.normalized(), v_des.normalized())) * 180.0/Pi;

                QPair<double, double> secondElevErrors = getSeconStepElevError(delta_deg, v_elev_deg, th1_deg, th2_deg);
                QPair<double, double> secondAzimuthErrors = getSeconStepAzimuthError(delta_deg, v_elev_deg, th1_deg, th2_deg);

                if(secondElevErrors.first != (double)-999.0 && secondAzimuthErrors.first != (double)-999.0){

                    QVector3D new_v_pol(v.x(), v.y() - (secondElevErrors.first * Pi/180.0), v.z() + (secondAzimuthErrors.first * Pi/180.0)); /*theta measured from vertical*/
                    QVector3D firstSecond_sun;
                    firstSecond_sun.setX(new_v_pol.x());

                    //                    for(double el = -secondElevErrors.second; el <= secondElevErrors.second; el+= 1.0){
                    //                        for(double az = -secondAzimuthErrors.second; az < secondAzimuthErrors.second; az+= 1.0){
                    firstSecond_sun.setY(new_v_pol.y());// + (el * Pi/180.0));
                    firstSecond_sun.setZ(new_v_pol.z());// + (az * Pi/180.0));

                    /*third step errors*/
                    double thirdError = 0;
                    double firstSecond_elev_deg = (Pi/2.0 - firstSecond_sun.y()) * 180.0/Pi;
                    if(bw == "best")
                        thirdError = 0.112482*pow(firstSecond_elev_deg, 0.784039);
                    if(bw == "second")
                        thirdError = 1.22733*pow(firstSecond_elev_deg, 0.353381);
                    if(bw == "worst")
                        thirdError = 0.0256134*pow(firstSecond_elev_deg, 1.58538);
                    if(thirdError < 0)
                        thirdError = 0;

                    QVector3D sunWithErrors;
                    sunWithErrors.setX(firstSecond_sun.x());
                    sunWithErrors.setZ(firstSecond_sun.z());

                    for(double elev = -thirdError; elev <= thirdError; elev += elev_res_deg){

                        sunWithErrors.setY(firstSecond_sun.y() + elev*Pi/180.0);
                        if((Pi/2.0 - sunWithErrors.y()) < 0) {
                            elev = -firstSecond_sun.y()*180.0/Pi;
                            sunWithErrors.setY(firstSecond_sun.y() + elev*Pi/180.0);
                        }


                        projected = QVector2D( tan(sunWithErrors.y())*cos(-sunWithErrors.z()) , tan(sunWithErrors.y())*sin(-sunWithErrors.z()) );

                        if(firstSecond_sun.y()<Pi/2.0-Pi/16  && projected.length()>R_MIN){ // horizonthoz ne legyen tul kozel (nagyon hosszu arnyek)

                            sShadows.append(projected);
                            paint = sShadows.last();
                            paint = transform.rotate2D(paint, North);
                            toPaint.append(paint);
                        }
                        QApplication::processEvents();
                    }
                }
                //                    }
                //                }
            }
        }
    }
//    emit signalWriteToList("Calculating errors ready.");
}

void AllSteps::calculateNorthErrors()
{
    QVector<double> nerrors;
    double NERROR;
    QVector<QVector2D> minimums;
    QVector2D INTERSECT;
    int h;

    int sShadowsSize = sShadows.size();
    int hypPointsSize = hypPoints.size();

    for(int i=0; i<sShadowsSize; i++){
        minimums.clear();
        nerrors.clear();

        QVector2D current_shadow = sShadows.at(i);
        sShadows[i] = transform.rotate2D(current_shadow, North);

        for (int j=1; j<hypPointsSize-1; j++){
            if( (fabs(sShadows.at(i).length()-hypPoints.at(j).length())<fabs(sShadows.at(i).length()-hypPoints.at(j-1).length()))&&(fabs(sShadows.at(i).length()-hypPoints.at(j).length())<fabs(sShadows.at(i).length()-hypPoints.at(j+1).length()))){
                minimums.append(hypPoints.at(j));
                nerrors.append(acos(QVector2D::dotProduct(hypPoints.at(j).normalized(), sShadows.at(i).normalized())));
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

            if (QVector3D::crossProduct(sShadows.at(i).toVector3D(), INTERSECT.toVector3D()).z() < 0)
                NERROR *= (-1.0);

            NERROR = NERROR/Pi*180.0;
            h = (int)(NERROR + 180) % 360;
            if(h >= 0){
                hist[h]++;
                hist[360]++;
            }
        }
    }
}

void AllSteps::paint()
{
    int size = 1355;
    QPainter painter;
    QPen pen;

    QImage im=QImage(size,size,QImage::Format_ARGB32);
    painter.begin(&im);

    painter.fillRect(0,0,size,size,Qt::white);
    painter.drawEllipse(QPointF(size/2.0,size/2.0),3,3);
    painter.drawEllipse(QPointF(size/2.0,size/2.0),5,5);
    painter.drawEllipse(QPointF(size/2.0,size/2.0),R_MIN*100, R_MIN*100);

    pen.setColor(Qt::red);
    painter.setPen(pen);
    int hypPointsSize = hypPoints.size();
    for(int q=0; q<hypPointsSize; q++)
        painter.drawPoint(size/2.0 + hypPoints.at(q).x()*100, size/2.0 - hypPoints.at(q).y()*100);

    pen.setColor(Qt::black);
    painter.setPen(pen);
    int toPaintSize = toPaint.size();
    for(int q=0; q<toPaintSize; q++)
        painter.drawPoint(size/2.0 + toPaint.at(q).x()*100, size/2.0 - toPaint.at(q).y()*100);

    painter.end();
    im.save("a.png");
}

void AllSteps::loadSecondErrorData()
{
    secondErrorElevList.clear();
    secondErrorAzimuthList.clear();
    paramRange.clear();

    QFile file("../result_ranges.csv");
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        qDebug("cannot open");
    QTextStream stream(&file);
    QString str;
    while(!stream.atEnd()){
        QString line = stream.readLine();
        if(!line.startsWith("#")){
            QStringList lineList = line.split("\t");
            str = lineList.takeFirst();
            secondErrorElevList[str].first = QString(lineList.takeFirst()).toDouble();
            secondErrorElevList[str].second = QString(lineList.takeFirst()).toDouble();
            secondErrorAzimuthList[str].first = QString(lineList.takeFirst()).toDouble();
            secondErrorAzimuthList[str].second = QString(lineList.takeFirst()).toDouble();
        }
    }

    for(int j = 0; j < secondErrorAzimuthList.keys().size(); j++){
        QString currentkey = secondErrorAzimuthList.keys().at(j);
        QStringList keyparts = currentkey.split(" ");

        for(int i = 0; i < keyparts.size(); i++) {
            QStringList currentPart = QString(keyparts.at(i)).split("-");
            double min = QString(currentPart.takeFirst()).toDouble();
            double max = QString(currentPart.takeFirst()).toDouble();
            if(i == 0){
                paramRange["delta_" + QString::number(j)].first = min;
                paramRange["delta_" + QString::number(j)].second = max;
            }
            if(i == 1){
                paramRange["elev_" + QString::number(j)].first = min;
                paramRange["elev_" + QString::number(j)].second = max;
            }
            if(i == 2){
                paramRange["th1_" + QString::number(j)].first = min;
                paramRange["th1_" + QString::number(j)].second = max;
            }
            if(i == 3){
                paramRange["th2_" + QString::number(j)].first = min;
                paramRange["th2_" + QString::number(j)].second = max;
            }
        }
    }

    secondErrorElevListKeysNum = secondErrorElevList.keys().size();
    secondErrorAzimuthListKeysNum = secondErrorAzimuthList.keys().size();

    QApplication::processEvents();
    emit signalWriteToList("Second error data loaded.");

}

QPair<double, double> AllSteps::getSeconStepAzimuthError(double delta, double elev, double th1, double th2)
{
    QPair<double, double> azimuthResult(-999, -999); /*to be easy to recognise*/
    for(int j = 0; j < secondErrorAzimuthListKeysNum; j++){
        QString currentkey = secondErrorAzimuthList.keys().at(j);
        if(delta >= paramRange["delta_" + QString::number(j)].first && delta <= paramRange["delta_" + QString::number(j)].second && elev >= paramRange["elev_" + QString::number(j)].first && elev <= paramRange["elev_" + QString::number(j)].second && th1 >= paramRange["th1_" + QString::number(j)].first && th1 <= paramRange["th1_" + QString::number(j)].second && th2 >= paramRange["th2_" + QString::number(j)].first && th2 <= paramRange["th2_" + QString::number(j)].second){
            azimuthResult.first = secondErrorAzimuthList[currentkey].first;
            azimuthResult.second = secondErrorAzimuthList[currentkey].second;
            break;
        }
    }

    return azimuthResult;
}

QPair<double, double> AllSteps::getSeconStepElevError(double delta, double elev, double th1, double th2)
{
    QPair<double, double> elevResult(-999, -999);
    for(int j = 0; j < secondErrorElevListKeysNum; j++){
        QString currentkey = secondErrorElevList.keys().at(j);
        if(delta >= paramRange["delta_" + QString::number(j)].first && delta <= paramRange["delta_" + QString::number(j)].second && elev >= paramRange["elev_" + QString::number(j)].first && elev <= paramRange["elev_" + QString::number(j)].second && th1 >= paramRange["th1_" + QString::number(j)].first && th1 <= paramRange["th1_" + QString::number(j)].second && th2 >= paramRange["th2_" + QString::number(j)].first && th2 <= paramRange["th2_" + QString::number(j)].second){
            elevResult.first = secondErrorElevList[currentkey].first;
            elevResult.second = secondErrorElevList[currentkey].second;
            break;
        }
    }

    return elevResult;
}


