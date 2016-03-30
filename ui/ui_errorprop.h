/********************************************************************************
** Form generated from reading UI file 'errorprop.ui'
**
** Created by: Qt User Interface Compiler version 5.6.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_ERRORPROP_H
#define UI_ERRORPROP_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>
#include <propagationform.h>

QT_BEGIN_NAMESPACE

class Ui_ErrorProp
{
public:
    QWidget *centralWidget;
    QGridLayout *gridLayout_2;
    PropagationForm *widget;
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *ErrorProp)
    {
        if (ErrorProp->objectName().isEmpty())
            ErrorProp->setObjectName(QStringLiteral("ErrorProp"));
        ErrorProp->resize(400, 300);
        centralWidget = new QWidget(ErrorProp);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        gridLayout_2 = new QGridLayout(centralWidget);
        gridLayout_2->setSpacing(6);
        gridLayout_2->setContentsMargins(11, 11, 11, 11);
        gridLayout_2->setObjectName(QStringLiteral("gridLayout_2"));
        widget = new PropagationForm(centralWidget);
        widget->setObjectName(QStringLiteral("widget"));

        gridLayout_2->addWidget(widget, 0, 0, 1, 1);

        ErrorProp->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(ErrorProp);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 400, 25));
        ErrorProp->setMenuBar(menuBar);
        mainToolBar = new QToolBar(ErrorProp);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        ErrorProp->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(ErrorProp);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        ErrorProp->setStatusBar(statusBar);

        retranslateUi(ErrorProp);

        QMetaObject::connectSlotsByName(ErrorProp);
    } // setupUi

    void retranslateUi(QMainWindow *ErrorProp)
    {
        ErrorProp->setWindowTitle(QApplication::translate("ErrorProp", "ErrorProp", 0));
    } // retranslateUi

};

namespace Ui {
    class ErrorProp: public Ui_ErrorProp {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_ERRORPROP_H
