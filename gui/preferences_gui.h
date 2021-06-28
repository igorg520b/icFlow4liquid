#ifndef PREFERENCES_GUI_H
#define PREFERENCES_GUI_H

#include <QObject>
#include <QString>
#include <QSettings>

// preferences related to visually presenting the data

class PreferencesGUI : public QObject
{
    Q_OBJECT
    Q_PROPERTY(bool LoadLastScn MEMBER LoadLastScene NOTIFY propertyChanged)

public:

    bool LoadLastScene = true;
    QString LastSceneFilename;  // filename for the last scene
    QString LastFolderFloes;    // folder from which last 2D flow was imported
    QString LastFolder3DGometry;    // floder from which last 3D object was imported

    int VisualizationOption = 0;
    bool ShowAxes = false;
    bool ShowScalarBar = true;
    bool ShowTentative = false;
    bool DrawEdges = true;
    bool DrawBoundary = true;
    bool DrawArrows = true;
    bool ShowWaterLevel = true;

    void SaveState(QSettings &settings)
    {
        settings.setValue("PreferencesGUI.LoadLastScene",LoadLastScene);
        settings.setValue("PreferencesGUI.LastSceneFilename",LastSceneFilename);
        settings.setValue("PreferencesGUI.LastFolderFloes",LastFolderFloes);
        settings.setValue("PreferencesGUI.LastFolder3DGometry",LastFolder3DGometry);

        settings.setValue("PreferencesGUI.VisualizationOption",VisualizationOption);
        settings.setValue("PreferencesGUI.ShowAxes",ShowAxes);
        settings.setValue("PreferencesGUI.ShowScalarBar",ShowScalarBar);
        settings.setValue("PreferencesGUI.ShowTentative",ShowTentative);
        settings.setValue("PreferencesGUI.DrawEdges",DrawEdges);
        settings.setValue("PreferencesGUI.DrawBoundary",DrawBoundary);
        settings.setValue("PreferencesGUI.DrawArrows",DrawArrows);
        settings.setValue("PreferencesGUI.ShowWaterLevel",ShowWaterLevel);
    }

    void LoadState(QSettings &settings)
    {
        LoadLastScene = settings.value("PreferencesGUI.LoadLastScene").toBool();
        LastSceneFilename = settings.value("PreferencesGUI.LastSceneFilename").toString();
        LastFolderFloes = settings.value("PreferencesGUI.LastFolderFloes").toString();
        LastFolder3DGometry = settings.value("PreferencesGUI.LastFolder3DGometry").toString();

        VisualizationOption = settings.value("PreferencesGUI.VisualizationOption").toInt();
        ShowAxes = settings.value("PreferencesGUI.ShowAxes").toBool();
        ShowScalarBar = settings.value("PreferencesGUI.ShowScalarBar").toBool();
        ShowTentative = settings.value("PreferencesGUI.ShowTentative").toBool();
        DrawEdges = settings.value("PreferencesGUI.DrawEdges").toBool();
        DrawBoundary = settings.value("PreferencesGUI.DrawBoundary").toBool();
        DrawArrows = settings.value("PreferencesGUI.DrawArrows").toBool();
        ShowWaterLevel = settings.value("PreferencesGUI.ShowWaterLevel").toBool();
    }

signals:
    void propertyChanged();
};

#endif // PREFERENCES_GUI_H
