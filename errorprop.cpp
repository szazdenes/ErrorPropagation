#include "errorprop.h"
#include "ui_errorprop.h"

ErrorProp::ErrorProp(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::ErrorProp)
{
    ui->setupUi(this);
}

ErrorProp::~ErrorProp()
{
    delete ui;
    QCoreApplication::quit();
}
