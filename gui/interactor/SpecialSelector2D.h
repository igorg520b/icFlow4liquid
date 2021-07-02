/**
 *
 * Camera rotation is not allowed with this interactor style.
 * Zooming affects the camera's parallel scale only, and assumes
 * that the camera is in parallel projection mode.
 * All camera changes invoke StartInteractionEvent when the button
 * is pressed, InteractionEvent when the mouse (or wheel) is moved,
 * and EndInteractionEvent when the button is released.  The bindings
 * are as follows:
 * Left mouse - Select/unselect a node.
 * Right mouse - Pin/unpin a node.
 * Middle mouse - Pan.
 * Scroll wheel - Zoom.
 */

#ifndef SPECIALSELECTOR2D_H
#define SPECIALSELECTOR2D_H

//#include "vtkInteractionStyleModule.h" // For export macro
#include <vtkInteractorStyle.h>
#include <vtkInteractorStyleRubberBand2D.h>
#include <vtkActor.h>
#include <vtkNew.h>
#include <vtkArcSource.h>
//#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>

class vtkUnsignedCharArray;

class MainWindow;

class SpecialSelector2D : public vtkInteractorStyleRubberBand2D
{
public:
    static SpecialSelector2D* New();
    vtkTypeMacro(SpecialSelector2D, vtkInteractorStyleRubberBand2D);

    SpecialSelector2D();

    MainWindow* mw;
    vtkNew<vtkArcSource> arcSource;
    vtkNew<vtkPolyDataMapper> mapper;
    vtkNew<vtkActor> actor; // selection circle

    void OnLeftButtonDown() override;
    void OnLeftButtonUp() override;
    void OnRightButtonDown() override;
    void OnRightButtonUp() override;
    void OnMouseMove() override;

    //  void OnMouseWheelForward() override;
    //  void OnMouseWheelBackward() override;

//  void OnMiddleButtonDown() override;
//  void OnMiddleButtonUp() override;

private:
    double selectionRadius = 0.25;
    double centerX = 0;
    double centerY = 0;

    void UpdateActor();


};

#endif
