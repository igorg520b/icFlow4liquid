#ifndef OBJECTPROPERTYBROWSER_H
#define OBJECTPROPERTYBROWSER_H

#include <QObject>
#include <QDebug>
#include <QPushButton>
#include <QMap>
#include <QMetaProperty>
#include "qteditorfactory.h"
#include "qttreepropertybrowser.h"
#include "qtpropertymanager.h"
#include "qtvariantproperty.h"

class ObjectPropertyBrowser : public QtTreePropertyBrowser
{
    Q_OBJECT

public:
    ObjectPropertyBrowser(QWidget* parent);
    void setActiveObject(QObject *obj);

private:
    QtVariantPropertyManager *variantManager;
    QObject *currentlyConnectedObject = nullptr;
    QMap<QtProperty *, const char*> propertyMap;


private slots:
    void valueChanged(QtProperty *property, const QVariant &value);

public slots:
    void objectUpdated();

};

#endif // OBJECTPROPERTYBROWSER_H
