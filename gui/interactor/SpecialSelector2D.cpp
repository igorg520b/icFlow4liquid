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

void SpecialSelector2D::OnLeftButtonDown()
{
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

void SpecialSelector2D::OnMouseMove()
{
    mouse_remained_stationary = false;
    //    std::cout << "SpecialSelector2D::OnMouseMove()" << std::endl;

    if(this->Interaction == SELECTING)
    {
        /*
        vtkRenderWindowInteractor* rwi = this->GetInteractor();
        int lastPt[] = { 0, 0 };
        rwi->GetLastEventPosition(lastPt);
        int curPt[] = { 0, 0 };
        rwi->GetEventPosition(curPt);

        if(this->CurrentRenderer == nullptr) std::cout << "CurrentRenderer is nullptr" << std::endl;

        vtkCamera* camera = this->CurrentRenderer->GetActiveCamera();

        vtkRenderer* renderer = this->CurrentRenderer;

        double camera_parallelScale = camera->GetParallelScale();

        int* renderer_getSize = renderer->GetSize();

        int renderer_getSize1 = renderer_getSize[1];

        double lastScale = 2.0 *  camera_parallelScale / renderer_getSize1;

        double lastFocalPt[] = { 0, 0, 0 };
        camera->GetFocalPoint(lastFocalPt);
        double lastPos[] = { 0, 0, 0 };
        camera->GetPosition(lastPos);

        double delta[] = { 0, 0, 0 };
        delta[0] = -lastScale * (curPt[0] - lastPt[0]);
        delta[1] = -lastScale * (curPt[1] - lastPt[1]);
        delta[2] = 0;

        for(icy::Node &nd : mw->modelController.model.mesh.nodes)
        {
            if(!nd.selected) continue;
            nd.intended_position.x() -= delta[0];
            nd.intended_position.y() -= delta[1];
            nd.vn = Eigen::Vector2d::Zero();
        }

        mw->modelController.model.UnsafeUpdateGeometry();

        rwi->Render();
        */
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

void SpecialSelector2D::OnRightButtonUp() { }

