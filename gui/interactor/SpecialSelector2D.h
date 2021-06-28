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

class vtkUnsignedCharArray;

class MainWindow;

class SpecialSelector2D : public vtkInteractorStyleRubberBand2D
{
public:
    static SpecialSelector2D* New();
    vtkTypeMacro(SpecialSelector2D, vtkInteractorStyleRubberBand2D);
    MainWindow* mw;

    void OnLeftButtonDown() override;
    void OnLeftButtonUp() override;
    void OnRightButtonDown() override;
    void OnRightButtonUp() override;
    void OnMouseMove() override;

/*
  void OnMiddleButtonDown() override;
  void OnMiddleButtonUp() override;
  void OnMouseWheelForward() override;
  void OnMouseWheelBackward() override;
*/
private:
    bool mouse_remained_stationary = false;

};

#endif
