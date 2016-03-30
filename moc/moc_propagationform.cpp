/****************************************************************************
** Meta object code from reading C++ file 'propagationform.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../propagationform.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'propagationform.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_PropagationForm_t {
    QByteArrayData data[10];
    char stringdata0[194];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_PropagationForm_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_PropagationForm_t qt_meta_stringdata_PropagationForm = {
    {
QT_MOC_LITERAL(0, 0, 15), // "PropagationForm"
QT_MOC_LITERAL(1, 16, 20), // "signalFirstStepStart"
QT_MOC_LITERAL(2, 37, 0), // ""
QT_MOC_LITERAL(3, 38, 20), // "signalThirdStepStart"
QT_MOC_LITERAL(4, 59, 19), // "signalAllStepsStart"
QT_MOC_LITERAL(5, 79, 15), // "slotWriteToList"
QT_MOC_LITERAL(6, 95, 6), // "string"
QT_MOC_LITERAL(7, 102, 30), // "on_firstStepPushButton_clicked"
QT_MOC_LITERAL(8, 133, 30), // "on_thirdStepPushButton_clicked"
QT_MOC_LITERAL(9, 164, 29) // "on_allStepsPushButton_clicked"

    },
    "PropagationForm\0signalFirstStepStart\0"
    "\0signalThirdStepStart\0signalAllStepsStart\0"
    "slotWriteToList\0string\0"
    "on_firstStepPushButton_clicked\0"
    "on_thirdStepPushButton_clicked\0"
    "on_allStepsPushButton_clicked"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_PropagationForm[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       7,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       3,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    0,   49,    2, 0x06 /* Public */,
       3,    0,   50,    2, 0x06 /* Public */,
       4,    0,   51,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       5,    1,   52,    2, 0x0a /* Public */,
       7,    0,   55,    2, 0x08 /* Private */,
       8,    0,   56,    2, 0x08 /* Private */,
       9,    0,   57,    2, 0x08 /* Private */,

 // signals: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

 // slots: parameters
    QMetaType::Void, QMetaType::QString,    6,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void PropagationForm::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        PropagationForm *_t = static_cast<PropagationForm *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->signalFirstStepStart(); break;
        case 1: _t->signalThirdStepStart(); break;
        case 2: _t->signalAllStepsStart(); break;
        case 3: _t->slotWriteToList((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 4: _t->on_firstStepPushButton_clicked(); break;
        case 5: _t->on_thirdStepPushButton_clicked(); break;
        case 6: _t->on_allStepsPushButton_clicked(); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (PropagationForm::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&PropagationForm::signalFirstStepStart)) {
                *result = 0;
                return;
            }
        }
        {
            typedef void (PropagationForm::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&PropagationForm::signalThirdStepStart)) {
                *result = 1;
                return;
            }
        }
        {
            typedef void (PropagationForm::*_t)();
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&PropagationForm::signalAllStepsStart)) {
                *result = 2;
                return;
            }
        }
    }
}

const QMetaObject PropagationForm::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_PropagationForm.data,
      qt_meta_data_PropagationForm,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *PropagationForm::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *PropagationForm::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_PropagationForm.stringdata0))
        return static_cast<void*>(const_cast< PropagationForm*>(this));
    return QWidget::qt_metacast(_clname);
}

int PropagationForm::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 7)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 7;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 7)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 7;
    }
    return _id;
}

// SIGNAL 0
void PropagationForm::signalFirstStepStart()
{
    QMetaObject::activate(this, &staticMetaObject, 0, Q_NULLPTR);
}

// SIGNAL 1
void PropagationForm::signalThirdStepStart()
{
    QMetaObject::activate(this, &staticMetaObject, 1, Q_NULLPTR);
}

// SIGNAL 2
void PropagationForm::signalAllStepsStart()
{
    QMetaObject::activate(this, &staticMetaObject, 2, Q_NULLPTR);
}
QT_END_MOC_NAMESPACE
