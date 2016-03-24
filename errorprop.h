#ifndef ERRORPROP_H
#define ERRORPROP_H

#include <QMainWindow>

namespace Ui {
class ErrorProp;
}

class ErrorProp : public QMainWindow
{
    Q_OBJECT

public:
    explicit ErrorProp(QWidget *parent = 0);
    ~ErrorProp();

private:
    Ui::ErrorProp *ui;
};

#endif // ERRORPROP_H
