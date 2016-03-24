#ifndef THIRDSTEP_H
#define THIRDSTEP_H

#define Pi 3.1415926536

#include <QObject>
#include <QWidget>
#include <QVector2D>
#include <QVector3D>
#include <QFile>
#include <QTextStream>
#include <math.h>
#include <QPainter>
#include <QPen>
#include <QtAlgorithms>
#include <QApplication>

#include "geomtransform.h"

class ThirdStep : public QObject
{
    Q_OBJECT

public:
    ThirdStep();

signals:
    void signalWriteToList(QString string);

public slots:
    void slotThirdStepStart();

private:
    QVector<QVector2D> hypPoints, sShadows;
    QVector3D s;
    QPair<QVector3D, QVector3D> sunWithErrors;
    double North;
    double R_MIN;
    double gnomonHeight;
    GeomTransform transform;
    QList<QVector2D> toPaint;

    double detectNorth(double sunElevDeg);
    void readHyperbolaFromFile(QString filename);
    double getMinRadius(QVector<QVector2D> hyperbola);
    void paint();
    void detectGnomonHeight(QString hypFile, double maxElevDeg);
    QVector2D projectSunShadow(double sunElevDeg);
    void calculateNErrors(QString date, QString ampm, QString bw);
};

#endif // THIRDSTEP_H
