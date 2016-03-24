#ifndef ALLSTEPS_H
#define ALLSTEPS_H

#define Pi 3.1415926536

#include <QObject>
#include <QWidget>

#include <QVector2D>
#include <QVector3D>
#include <QFileDialog>
#include <QDir>
#include <QTextStream>
#include <math.h>
#include <QPainter>
#include <QPen>
#include <QtAlgorithms>
#include <QApplication>
#include <QImage>
#include <QColor>
#include <QPainter>
#include <QTime>
#include <QFuture>
#include <QtConcurrent/QtConcurrent>

#include "geomtransform.h"

class AllSteps : public QObject
{
    Q_OBJECT
public:
    explicit AllSteps(QObject *parent = 0);

signals:
    void signalWriteToList(QString string);

public slots:
    void slotAllStepsStart();

private:
    QVector<QVector2D> hypPoints, sShadows, toPaint;
    QVector<QVector3D> p2s;
    QVector<double> deg_p2s;
    QVector<int> hist;
    QVector3D p1, p2, s;

    double North;
    double R_MIN;
    double deg_p1, deg_p2, err1, err2;
    QImage *image;
    GeomTransform transform;
//    QList<QVector2D> toPaint;
    int size;

    /*from thirdstep*/
//    double detectNorth(double sunElevDeg);
    double getMinRadius(QVector<QVector2D> hyperbola);
//    void calculateNErrors(QString date, QString ampm, QString bw);

    /*from firststep*/
    void initializeHist();
    double getError(int mode, double deg_of_pol, int stein);
    void readHyperbolaFromFile(QString filename);
    void writeDatFile(QString imname, QString hyperbola, QString timeOfDay, int stein, int p1_res, double p2_res, double deltoid_res, int mode, QString bestworst);
    void readImage(QString filename);
    void detectSun();
    void detectNorth();
    void calculatep2s(QVector3D p1, double min, double max, double res);
    QVector3D intersectionOfGreatCircles(QVector3D GC1A, QVector3D GC1B, QVector3D GC2A, QVector3D GC2B);
    void calculateSunPositions(QVector3D sun, QVector3D point1, QVector3D point2, double error1, double error2, double res, double elev_res_deg, QString bw);
    void calculateNorthErrors();
    void paint();

    /*second step*/
    void loadSecondErrorData();
    QPair<double, double> getSeconStepElevError(double delta, double elev, double th1, double th2); //first:ave, second:std
    QPair<double, double> getSeconStepAzimuthError(double delta, double elev, double th1, double th2);
    QMap<QString, QPair<double, double> > secondErrorElevList, secondErrorAzimuthList;
    QMap<QString, QPair<double, double> > paramRange;
};

#endif // ALLSTEPS_H
