#include "thirdstep.h"

ThirdStep::ThirdStep()
{

}

void ThirdStep::slotThirdStepStart()
{
    calculateNErrors("0621", "am", "best");
    calculateNErrors("0621", "am", "secondworst");
    calculateNErrors("0621", "am", "worst");
    calculateNErrors("0621", "pm", "best");
    calculateNErrors("0621", "pm", "secondworst");
    calculateNErrors("0621", "pm", "worst");
    calculateNErrors("0922", "am", "best");
    calculateNErrors("0922", "am", "secondworst");
    calculateNErrors("0922", "am", "worst");
    calculateNErrors("0922", "pm", "best");
    calculateNErrors("0922", "pm", "secondworst");
    calculateNErrors("0922", "pm", "worst");
//    paint();
}

double ThirdStep::detectNorth(double sunElevDeg)
{
    double NError;
    QVector<double> nerrors;
    QVector<QVector2D> minimums;
    QVector2D INTERSECT;

    QVector2D sun_shadow = projectSunShadow(sunElevDeg);

    minimums.clear();
    nerrors.clear();

    int hypPointsSize = hypPoints.size();

    if(sun_shadow.length() >= R_MIN){
        for (int j=1; j<hypPointsSize-1; j++){
            if( (fabs(sun_shadow.length()-hypPoints.at(j).length())<fabs(sun_shadow.length()-hypPoints.at(j-1).length()))&&(fabs(sun_shadow.length()-hypPoints.at(j).length())<fabs(sun_shadow.length()-hypPoints.at(j+1).length()))){
                minimums.append(hypPoints.at(j));
                nerrors.append(acos(QVector2D::dotProduct(hypPoints.at(j).normalized(), sun_shadow.normalized())));
            }
        }
    }

    if(minimums.isEmpty() || nerrors.isEmpty()){
        emit signalWriteToList("Error: Too short sun shadow length.");
        return 999;
    }

    INTERSECT = minimums.at(0);
    NError = nerrors.at(0);

    int minimumsSize = minimums.size();
    for(int j=0; j<minimumsSize; j++)
        if(nerrors.at(j) < NError){
            INTERSECT = minimums.at(j);
            NError = nerrors.at(j);
        }


    if (QVector3D::crossProduct(QVector3D (sun_shadow.x(), sun_shadow.y(),0.0), QVector3D (INTERSECT.x(), INTERSECT.y(), 0.0) ).z() < 0)
        NError *= (-1.0);
    return NError;
}

void ThirdStep::readHyperbolaFromFile(QString filename)
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

double ThirdStep::getMinRadius(QVector<QVector2D> hyperbola)
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

void ThirdStep::paint()
{
    QPainter painter;
    QPen pen;
    int size = 1000;

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


    pen.setWidth(5);
    pen.setColor(Qt::black);
    int toPaintSize = toPaint.size();
    for(int q=0; q<toPaintSize; q++){
        if(q==1)
            pen.setColor(Qt::green);
        if(q==2)
            pen.setColor(Qt::blue);
        if(q==4)
            pen.setColor(Qt::green);
        if(q==5)
            pen.setColor(Qt::blue);
        painter.setPen(pen);
        painter.drawPoint(size/2.0 + toPaint.at(q).x()*100, size/2.0 - toPaint.at(q).y()*100);
    }

    painter.end();
    im.save("a.png");
}

void ThirdStep::detectGnomonHeight(QString hypFile, double maxElevDeg)
{
    QList<QVector2D> hyperbola;

    QFile file(hypFile);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
        qDebug("baj");

    QTextStream filestream(&file);
    double x,y;

    while (!filestream.atEnd()){
        QString line = filestream.readLine();
        QTextStream stream(&line);
        stream >> x >> y;
        hyperbola.append(QVector2D(x,y));
    }

    QList<double> hypLength;

    for (int i = 0; i < hyperbola.size(); i++)
        hypLength.append(hyperbola.at(i).length());


    file.close();

    qSort(hypLength.begin(), hypLength.end());
    R_MIN = hypLength.first();
    gnomonHeight = R_MIN*tan(maxElevDeg*Pi/180.0);

    emit signalWriteToList("Gnomon height: " + QString::number(gnomonHeight));
}

QVector2D ThirdStep::projectSunShadow(double sunElevDeg)
{
    QVector2D shadowPos;
    shadowPos.setX(gnomonHeight/tan(sunElevDeg*Pi/180.0));
    shadowPos.setY(0);
    return shadowPos;
}

void ThirdStep::calculateNErrors(QString date, QString ampm, QString bw)
{
    double maxElev;
    QString hypname;
    if(date == "0621"){
        hypname = "/home/denes/Documents/Labor/Viking/ErrorPropagation/l61d0621.dat";
        maxElev = 52.5;
    }
    if(date == "0922"){
        hypname = "/home/denes/Documents/Labor/Viking/ErrorPropagation/l61d0922.dat";
        maxElev = 29;
    }
    detectGnomonHeight(hypname, maxElev);
    hypname.remove(".dat").append(ampm + ".dat");
    readHyperbolaFromFile(hypname);

    double error = 0, NE1, NE2, estNorth1, estNorth2;
    QFile outfile("/home/denes/Documents/Labor/Viking/ErrorPropagation/step3results/step3Error" + QString("_") + date + QString("_") + ampm + QString("_") + bw + ".csv");
    if (!outfile.open(QIODevice::WriteOnly | QIODevice::Text))
        qDebug("baj");

    QTextStream out(&outfile);
    out << "elevation\terror1\terror2\n";

    for(double i = 0.1; i < maxElev; i+=0.1){
        if(bw == "best")
            error = 0.112482*pow(i, 0.784039);
        if(bw == "secondworst")
            error = 1.22733*pow(i, 0.353381);
        if(bw == "worst")
            error = 0.0256134*pow(i, 1.58538);
        if(error < 0)
            error = 0;
        North = detectNorth(i);
        emit signalWriteToList("north detected: " + QString::number(North*180/Pi) + " deg");
        if(i-error<0)
            estNorth1 = detectNorth(i);
        else
            estNorth1 = detectNorth(i-error);
        estNorth2 = detectNorth(i+error);


        NE1 = North*180/Pi - estNorth1*180/Pi;
        NE2 = North*180/Pi - estNorth2*180/Pi;

        if(qAbs(NE1) > 360)
            out << QString::number(i) + "\t" + "NaN" + "\t" + QString::number(NE2) + "\n";
        else if(qAbs(NE2) > 360)
            out << QString::number(i) + "\t" + QString::number(NE1) + "\t" + "NaN" + "\n";
        else if(qAbs(NE1) > 360 && qAbs(NE2) > 360)
            out << QString::number(i) + "\t" + QString::number(NE1) + "\t" + "NaN" + "\n";
        else
            out << QString::number(i) + "\t" + QString::number(NE1) + "\t" + QString::number(NE2) + "\n";
        QApplication::processEvents();
    }

    emit signalWriteToList(outfile.fileName() + " ready");
    outfile.close();

}
