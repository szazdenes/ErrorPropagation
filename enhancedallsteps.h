#ifndef ENHANCEDALLSTEPS_H
#define ENHANCEDALLSTEPS_H

#include <QObject>
#include <QWidget>

class EnhancedAllsteps : public QObject
{
    Q_OBJECT
public:
    explicit EnhancedAllsteps(QObject *parent = 0);

signals:
    void signalWriteToList(QString string);

public slots:
    void slotEnhAllStepsStart();
};

#endif // ENHANCEDALLSTEPS_H
