#ifndef ENHANCEDALLSTEPS_H
#define ENHANCEDALLSTEPS_H

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

class EnhancedAllsteps : public QObject
{
    Q_OBJECT
public:
    explicit EnhancedAllsteps(QObject *parent = 0);

signals:
    void signalWriteToList(QString string);

public slots:
    void slotEnhAllStepsStart();

private:
    GeomTransform transform;

    void paint(QList<QVector2D> hyp_points, QList<QVector2D> to_paint, double r_min);
    QImage readImage(QString filename);

    double getMinRadius(QList<QVector2D> hyperbola);
    QList<QVector2D> readHyperbolaFromFile(QString filename);
    double getError(double deg_of_pol, int stein);
    void writeDatFile(QString imname, QString hyperbola, QString timeOfDay, int stein, int p1_res, double p2_res, double deltoid_res, double second_error_res, double elev_res, QVector<int> histogram, double north);
    QVector3D detectSun(QImage im);
    double detectNorth(QList<QVector2D> hyp_points, QVector3D sun);
    QPair<QList<QVector3D>, QList<double> > calculatep2s(QImage im, QVector3D sun, QVector3D p1, double min, double max, double res);
    QVector3D intersectionOfGreatCircles(QVector3D GC1A, QVector3D GC1B, QVector3D GC2A, QVector3D GC2B, QVector3D sun);
    QVector<int> calculateNorthErrors(QList<QVector2D> sun_shadows, QList<QVector2D> hyp_points, double north);

    QList<QPair<QPair<QVector3D, QVector3D>, QVector3D> > getIntersectionsFromFirstStep_pol(QList<QVector3D> sunAndPoints, QList<double> firstStepErrors, QVector3D sun); //sunAndPoint 0: sun 1: point1 2: point2, firstStepErrors 0: error1 1: error2 2: deltoid_res
    QList<QVector3D> getPointsFromSecondStep_pol(QList<QPair<QPair<QVector3D, QVector3D>, QVector3D> > *intersectPoints_pol, double seconderror_resolution);
    QList<QVector3D> getPointsFromThirdStep_pol(QList<QVector3D> *secondErrorPoint_pol, double elev_res_deg);
    QPair<QList<QVector2D>, QList<QVector2D> > getSunShadows(QList<QVector3D> *thirdErrorPoints_pol, double r_min, double north); //first: sunShadow, second: toPaint

    /*second step*/
    void loadSecondErrorData();
    QPair<double, double> getSeconStepError(QList<double> params, QMap<QString, QPair<double, double> > secondErrorMap, QMap<QString, QPair<double, double> > rangeMap, double keysNum); //first:ave, second:std
    QMap<QString, QPair<double, double> > secondErrorElevList, secondErrorAzimuthList;
    QMap<QString, QPair<double, double> > paramRange;
    int secondErrorElevListKeysNum, secondErrorAzimuthListKeysNum;

};

#endif // ENHANCEDALLSTEPS_H
