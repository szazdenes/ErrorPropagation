/********************************************************************************
** Form generated from reading UI file 'propagationform.ui'
**
** Created by: Qt User Interface Compiler version 5.6.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_PROPAGATIONFORM_H
#define UI_PROPAGATIONFORM_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QListWidget>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_PropagationForm
{
public:
    QGridLayout *gridLayout_2;
    QPushButton *thirdStepPushButton;
    QPushButton *firstStepPushButton;
    QPushButton *secondStepPushButton;
    QPushButton *allStepsPushButton;
    QListWidget *listWidget;

    void setupUi(QWidget *PropagationForm)
    {
        if (PropagationForm->objectName().isEmpty())
            PropagationForm->setObjectName(QStringLiteral("PropagationForm"));
        PropagationForm->resize(703, 563);
        gridLayout_2 = new QGridLayout(PropagationForm);
        gridLayout_2->setObjectName(QStringLiteral("gridLayout_2"));
        thirdStepPushButton = new QPushButton(PropagationForm);
        thirdStepPushButton->setObjectName(QStringLiteral("thirdStepPushButton"));

        gridLayout_2->addWidget(thirdStepPushButton, 0, 2, 1, 1);

        firstStepPushButton = new QPushButton(PropagationForm);
        firstStepPushButton->setObjectName(QStringLiteral("firstStepPushButton"));

        gridLayout_2->addWidget(firstStepPushButton, 0, 0, 1, 1);

        secondStepPushButton = new QPushButton(PropagationForm);
        secondStepPushButton->setObjectName(QStringLiteral("secondStepPushButton"));

        gridLayout_2->addWidget(secondStepPushButton, 0, 1, 1, 1);

        allStepsPushButton = new QPushButton(PropagationForm);
        allStepsPushButton->setObjectName(QStringLiteral("allStepsPushButton"));

        gridLayout_2->addWidget(allStepsPushButton, 0, 3, 1, 1);

        listWidget = new QListWidget(PropagationForm);
        listWidget->setObjectName(QStringLiteral("listWidget"));

        gridLayout_2->addWidget(listWidget, 1, 0, 1, 4);


        retranslateUi(PropagationForm);

        QMetaObject::connectSlotsByName(PropagationForm);
    } // setupUi

    void retranslateUi(QWidget *PropagationForm)
    {
        PropagationForm->setWindowTitle(QApplication::translate("PropagationForm", "Form", 0));
        thirdStepPushButton->setText(QApplication::translate("PropagationForm", "Third step", 0));
        firstStepPushButton->setText(QApplication::translate("PropagationForm", "First step", 0));
        secondStepPushButton->setText(QApplication::translate("PropagationForm", "Second step", 0));
        allStepsPushButton->setText(QApplication::translate("PropagationForm", "All steps", 0));
    } // retranslateUi

};

namespace Ui {
    class PropagationForm: public Ui_PropagationForm {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_PROPAGATIONFORM_H
