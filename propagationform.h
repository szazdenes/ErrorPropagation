#ifndef PROPAGATIONFORM_H
#define PROPAGATIONFORM_H

#include <QWidget>

#include "firststep.h"
#include "thirdstep.h"
#include "allsteps.h"

namespace Ui {
class PropagationForm;
}

class PropagationForm : public QWidget
{
    Q_OBJECT

public:
    explicit PropagationForm(QWidget *parent = 0);
    ~PropagationForm();

signals:
    void signalFirstStepStart();
    void signalThirdStepStart();
    void signalAllStepsStart();

public slots:
    void slotWriteToList(QString string);

private slots:
    void on_firstStepPushButton_clicked();

    void on_thirdStepPushButton_clicked();

    void on_allStepsPushButton_clicked();

private:
    Ui::PropagationForm *ui;
};

#endif // PROPAGATIONFORM_H
