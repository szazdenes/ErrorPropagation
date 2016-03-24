#ifndef FIRSTSTEP_H
#define FIRSTSTEP_H

#include <QObject>
#include <QWidget>
#include <QImage>
#include <QColor>
#include <QPainter>
#include <QDebug>
#include <QVector2D>
#include <QVector3D>
#include <QFile>
#include <math.h>
#include <QCoreApplication>

#include <geomtransform.h>

#define Pi 3.1415926536

class FirstStep : public QObject
{
    Q_OBJECT

public:
    FirstStep();

signals:
    void signalWriteToList(QString string);

public slots:
    void slotFirstStepStart();

private:
    QVector<QVector2D> hypPoints, sShadows, toPaint;
    QImage *image;
    QVector<QVector3D> p2s;
    QVector<double> deg_p2s;
    QVector<int> hist;
    QVector3D p1, p2, s;
    double deg_p1, deg_p2, err1, err2;
    double North;
    double R_MIN;
    int size;
    GeomTransform transform;

    void initializeHist();
    double getError(int mode, double deg_of_pol, int stein);
    void readHyperbolaFromFile(QString filename);
    double getMinRadius(QVector<QVector2D> hyperbola);
    void writeDatFile(int no, QString hyperbola, QString timeOfDay, int stein, int p1_res, double p2_res, double deltoid_res, int mode);
    void readImage(QString filename);
    void detectSun();
    void detectNorth();
    void calculatep2s(QVector3D p1, double min, double max, double res);
    QVector3D intersectionOfGreatCircles(QVector3D GC1A, QVector3D GC1B, QVector3D GC2A, QVector3D GC2B);
    void calculateSunPositions(QVector3D sun, QVector3D point1, QVector3D point2, double error1, double error2, double res);
    void calculateNorthErrors();
    void paint();
};

#endif // FIRSTSTEP_H
