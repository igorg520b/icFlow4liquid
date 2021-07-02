#include "SpecialSelector2D.h"
#include "mainwindow.h"

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCommand.h>
#include <vtkObjectFactory.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPointPicker.h>
#include <vtkRendererCollection.h>

#include <iostream>

vtkStandardNewMacro(SpecialSelector2D);

SpecialSelector2D::SpecialSelector2D()
{
    this->Interaction = SELECTING;

    arcSource->SetResolution(100);
    arcSource->SetCenter(0,0,0);
    arcSource->NegativeOff();
    arcSource->UseNormalAndAngleOn();
    arcSource->SetNormal(0,0,1);
    arcSource->SetPolarVector(1,0,0);
    arcSource->SetAngle(360.0);


    arcSource->Update();
    mapper->SetInputConnection(arcSource->GetOutputPort());
//    mapper->SetInputConnection(sphereSource->GetOutputPort());
    actor->SetMapper(mapper);
    actor->GetProperty()->EdgeVisibilityOn();
    actor->GetProperty()->SetColor(0.1, 0.1, 0.1);
    actor->GetProperty()->SetEdgeColor(90.0/255.0, 90.0/255.0, 97.0/255.0);
    actor->GetProperty()->LightingOff();
    actor->GetProperty()->ShadingOff();
    actor->GetProperty()->SetInterpolationToFlat();
    actor->PickableOff();
    actor->GetProperty()->SetLineWidth(3);

    actor->SetVisibility(this->Interaction == SELECTING);
}

void SpecialSelector2D::GetCurrentCoords(double &x, double &y)
{
    vtkRenderWindowInteractor* rwi = this->GetInteractor();
    int curPt[] = { 0, 0 };
    rwi->GetEventPosition(curPt);

    vtkRenderer* renderer = rwi->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
    vtkCamera* camera = renderer->GetActiveCamera();

    double camera_parallelScale = camera->GetParallelScale();
    int* renderer_getSize = renderer->GetSize();
    int renderer_getSize1 = renderer_getSize[1];

    double camPosX, camPosY, camPosZ;
    camera->GetPosition(camPosX,camPosY,camPosZ);

    double lastScale = 2.0 *  camera_parallelScale / renderer_getSize1;

    x = lastScale * (curPt[0]-renderer_getSize[0]/2)+camPosX;
    y = lastScale * (curPt[1]-renderer_getSize[1]/2)+camPosY;
}


void SpecialSelector2D::UpdateActor()
{
    GetCurrentCoords(centerX, centerY);
    arcSource->SetCenter(centerX,centerY,0);
    arcSource->SetPolarVector(selectionRadius,0,0);
    arcSource->Update();
    this->GetInteractor()->Render();
}


void SpecialSelector2D::OnLeftButtonDown()
{
    if(this->Interaction == SELECTING)
    {
        stretching = true;
        GetCurrentCoords(initialX, initialY);
        mw->model.AttachSpring(initialX, initialY, selectionRadius);
    }
    /*
    this->FindPokedRenderer(this->StartPosition[0], this->StartPosition[1]);

    if(this->Interaction != NONE) return;
    this->Interaction = SELECTING; // this stands for selecting or moving
    // record mouse position
    this->StartPosition[0] = this->Interactor->GetEventPosition()[0];
    this->StartPosition[1] = this->Interactor->GetEventPosition()[1];

    // set "mouse_remained_statinary" flag
    mouse_remained_stationary = true;

    //    std::cout << "SpecialSelector2D::OnLeftButtonDown()" << std::endl;
    */
}

void SpecialSelector2D::OnLeftButtonUp()
{
    if(this->Interaction == SELECTING)
    {
        stretching = false;
        mw->model.ReleaseSpring();
    }

    /*
    if(this->Interaction != SELECTING) return;
    this->Interaction = NONE;
    // if "mouse_remained_statinary", then toggle point selection
    if(mouse_remained_stationary)
    {
        //        std::cout << "SpecialSelector2D::OnLeftButtonUp() mouse_remained_stationary" << std::endl;
        vtkPointPicker *pp = dynamic_cast<vtkPointPicker*>(this->Interactor->GetPicker());

        int result = pp->Pick(this->Interactor->GetEventPosition()[0],
                this->Interactor->GetEventPosition()[1],
                0,
                this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());

        vtkIdType id = pp->GetPointId();
        if(result == 0 || id<0) {
            //            std::cout << "SpecialSelector2D::OnLeftButtonUp() nothing selected" << std::endl;
            return;
        }

        mw->modelController.model.mesh.nodes[id].selected = !mw->modelController.model.mesh.nodes[id].selected;
        if(mw->modelController.model.mesh.nodes[id].selected)
        {
            mw->modelController.model.mesh.nodes[id].pinned = true;
            mw->modelController.model.mesh.nodes[id].intended_position = mw->modelController.model.mesh.nodes[id].xn;
        }
        mw->modelController.model.UnsafeUpdateGeometry();
        this->GetInteractor()->Render();
    }
    */
}


void SpecialSelector2D::OnMouseWheelForward()
{
    if(this->Interaction == SELECTING)
    {
        selectionRadius/=1.1;
        UpdateActor();
    }
    else
    {
        vtkInteractorStyleRubberBand2D::OnMouseWheelForward();
    }
}

void SpecialSelector2D::OnMouseWheelBackward()
{
    if(this->Interaction == SELECTING)
    {
        selectionRadius*=1.1;
        UpdateActor();
    }
    else
    {
        vtkInteractorStyleRubberBand2D::OnMouseWheelBackward();
    }
}



void SpecialSelector2D::OnMouseMove()
{
//    mouse_remained_stationary = false;
    //    std::cout << "SpecialSelector2D::OnMouseMove()" << std::endl;

    if(this->Interaction == SELECTING)
    {
        UpdateActor();
        if(stretching) mw->model.AdjustSpring(centerX-initialX, centerY-initialY);
    }
    else
    {
        vtkInteractorStyleRubberBand2D::OnMouseMove();
    }
}


void SpecialSelector2D::OnRightButtonDown()
{
    /*
    //vtkPointPicker
    vtkPointPicker *pp = dynamic_cast<vtkPointPicker*>(this->Interactor->GetPicker());

    int result = pp->Pick(this->Interactor->GetEventPosition()[0],
            this->Interactor->GetEventPosition()[1],
            0,
            this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer());

    vtkIdType id = pp->GetPointId();
    if(result == 0 || id<0) return;

    mw->modelController.model.mesh.nodes[id].pinned = !mw->modelController.model.mesh.nodes[id].pinned;
    if(mw->modelController.model.mesh.nodes[id].pinned == false)
        mw->modelController.model.mesh.nodes[id].selected = false;
    mw->modelController.model.UnsafeUpdateGeometry();
    this->GetInteractor()->Render();
    */
}

void SpecialSelector2D::OnRightButtonUp()
{
    if(this->Interaction == NONE) this->Interaction = SELECTING;
    else this->Interaction = NONE;
    actor->SetVisibility(this->Interaction == SELECTING);
    if(this->Interaction == SELECTING) UpdateActor();
    else this->GetInteractor()->Render();
}

