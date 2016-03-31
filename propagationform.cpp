#include "propagationform.h"
#include "ui_propagationform.h"

PropagationForm::PropagationForm(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::PropagationForm)
{
    ui->setupUi(this);
}

PropagationForm::~PropagationForm()
{
    delete ui;
}

void PropagationForm::slotWriteToList(QString string)
{
    if(ui->listWidget->count() > 100)
        ui->listWidget->takeItem(0);

    ui->listWidget->addItem(string);
    ui->listWidget->scrollToBottom();
}

void PropagationForm::on_firstStepPushButton_clicked()
{
    ui->listWidget->clear();
    FirstStep first;
    connect(this, &PropagationForm::signalFirstStepStart, &first, &FirstStep::slotFirstStepStart);
    connect(&first, &FirstStep::signalWriteToList, this, &PropagationForm::slotWriteToList);

    emit signalFirstStepStart();
}

void PropagationForm::on_thirdStepPushButton_clicked()
{
    ui->listWidget->clear();
    ThirdStep third;
    connect(this, &PropagationForm::signalThirdStepStart, &third, &ThirdStep::slotThirdStepStart);
    connect(&third, &ThirdStep::signalWriteToList, this, &PropagationForm::slotWriteToList);

    emit signalThirdStepStart();
}

void PropagationForm::on_allStepsPushButton_clicked()
{
    ui->listWidget->clear();
    AllSteps all;
    connect(this, &PropagationForm::signalAllStepsStart, &all, &AllSteps::slotAllStepsStart);
    connect(&all, &AllSteps::signalWriteToList, this, &PropagationForm::slotWriteToList);

    emit signalAllStepsStart();
}

void PropagationForm::on_pushButton_clicked()
{
    ui->listWidget->clear();
    EnhancedAllsteps enhall;
    connect(this, &PropagationForm::signalEnhancedAllStepsStart, &enhall, &EnhancedAllsteps::slotEnhAllStepsStart);
    connect(&enhall, &EnhancedAllsteps::signalWriteToList, this, &PropagationForm::slotWriteToList);

    emit signalEnhancedAllStepsStart();
}
