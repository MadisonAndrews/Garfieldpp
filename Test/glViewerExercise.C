//script showing how to use the GL viewer API to animate a picture
//Author: Richard maunder

#include "TGLViewer.h"
#include "TGLPerspectiveCamera.h"
#include "TTimer.h"
#include "TRandom.h"
#include "TVirtualPad.h"

TGLViewer::ECameraType camera;
TTimer timer(25);
TRandom randGen(0);

Int_t moveCount = 0;

void AnimateCamera()
{
   // initialization
   static Double_t fov = 30;
   static Double_t zoom = 0.78;
   static Double_t dolly = 1500.0;
   static Double_t center[3] = {-164.0, -164.0, -180.0};
   static Double_t hRotate = 0.0;
   static Double_t vRotate = 0.0;
   // steps
   static Double_t fovStep = randGen.Rndm()*3.0 - 0.5;
   static Double_t zoomStep = (20 - randGen.Rndm())/1000.;
   static Double_t dollyStep = randGen.Rndm()*5.0 - 1.0;
   static Double_t centerStep[3] = {randGen.Rndm()*4, randGen.Rndm()*4,
                                    randGen.Rndm()*4 };
   static Double_t hRotateStep = randGen.Rndm()*0.025;
   static Double_t vRotateStep = randGen.Rndm()*0.05;

   // move center
   center[0] += centerStep[0];
   center[1] += centerStep[1];
   center[2] += centerStep[2];
   Double_t mag = TMath::Sqrt(center[0]*center[0] + center[1]*center[1] +
                              center[2]*center[2]);
   if(mag > 500)
   {
      centerStep[0] = -centerStep[0];
      centerStep[1] = -centerStep[1];
      centerStep[2] = -centerStep[2];
   }

   // rotate
   hRotate += hRotateStep;
   vRotate += vRotateStep;
   if (vRotate >= TMath::TwoPi() || vRotate <= 0.0)
      vRotateStep = -vRotateStep;

   if (hRotate >= (TMath::PiOver2()- 0.02f) ||
       hRotate <= (0.07f -TMath::PiOver2())) {
      hRotateStep = -hRotateStep;
   }

   // dolly
   dolly += dollyStep;
   if (dolly >= 2000.0 || dolly <= 1500.0)
      dollyStep = -dollyStep;

   // modify frustum
   TGLViewer * v = (TGLViewer *)gPad->GetViewer3D();
   if(camera < 3)
   {
      fov += fovStep;
      if (fov > 130.0 || fov < 5.0)
         fovStep = - fovStep; }
   else
   {
      zoom += zoomStep;
      if (zoom > 4.0 || zoom < 0.25)
         zoomStep = - zoomStep;
   }

   // apply
   if(camera < 3)
      v->SetPerspectiveCamera(camera, fov, dollyStep, center, hRotateStep,
                              vRotateStep);
   else
      v->SetOrthoCamera(camera, zoom, dollyStep, center, hRotateStep,
                        vRotateStep);

   if (++moveCount % 10 == 0)
      v->RefreshPadEditor(v);
}

void glViewerExercise()
{
  //   gROOT->ProcessLine(".x nucleus.C");
   TGLViewer * v = (TGLViewer *)gPad->GetViewer3D();

   Double_t NeutronRadius = 60,
            ProtonRadius = 60,
            NucleusRadius,
            distance = 60;
   Double_t vol = nProtons + nNeutrons;
   vol = 3 * vol / (4 * TMath::Pi());

   NucleusRadius = distance * TMath::Power(vol, 1./3.);
//   cout << "NucleusRadius: " << NucleusRadius << endl;

   TGeoManager * geom = new TGeoManager("nucleus", "Model of a nucleus");
   geom->SetNsegments(40);
   TGeoMaterial *matEmptySpace = new TGeoMaterial("EmptySpace", 0, 0, 0);
   TGeoMaterial *matProton     = new TGeoMaterial("Proton"    , .938, 1., 10000.);
   TGeoMaterial *matNeutron    = new TGeoMaterial("Neutron"   , .935, 0., 10000.);

   TGeoMedium *EmptySpace = new TGeoMedium("Empty", 1, matEmptySpace);
   TGeoMedium *Proton     = new TGeoMedium("Proton", 1, matProton);
   TGeoMedium *Neutron    = new TGeoMedium("Neutron",1, matNeutron);

//  the space where the nucleus lives (top container volume)

   Double_t worldx = 200.;
   Double_t worldy = 200.;
   Double_t worldz = 200.;
 
   TGeoVolume *top = geom->MakeBox("WORLD", EmptySpace, worldx, worldy, worldz); 
   geom->SetTopVolume(top);

   TGeoVolume * proton  = geom->MakeSphere("proton",  Proton,  0., ProtonRadius); 
   TGeoVolume * neutron = geom->MakeSphere("neutron", Neutron, 0., NeutronRadius); 
   proton->SetLineColor(kRed);
   neutron->SetLineColor(kBlue);

   Double_t x, y, z, dummy;
   Int_t i = 0; 
   // animate : add random nodes

   /*
   // Random draw style
   Int_t style = randGen.Integer(3);
   switch (style)
   {
      case 0: v->SetStyle(TGLRnrCtx::kFill); break;
      case 1: v->SetStyle(TGLRnrCtx::kOutline); break;
      case 2: v->SetStyle(TGLRnrCtx::kWireFrame); break;
   }

   // Lights - turn some off randomly
   TGLLightSet* ls = v->GetLightSet();
   if (randGen.Integer(2) == 0)
      ls->SetLight(TGLLightSet::kLightLeft, kFALSE);
   if (randGen.Integer(2) == 0)
      ls->SetLight(TGLLightSet::kLightRight, kFALSE);
   if (randGen.Integer(2) == 0)
      ls->SetLight(TGLLightSet::kLightTop, kFALSE);
   if (randGen.Integer(2) == 0)
      ls->SetLight(TGLLightSet::kLightBottom, kFALSE);

   // Random camera type
   Int_t id = randGen.Integer(6);
   camera = (TGLViewer::ECameraType)id;
   v->SetCurrentCamera(camera);
   v->CurrentCamera().SetExternalCenter(kTRUE);
   if (id > 2) {
      TGLOrthoCamera& o = v->CurrentCamera();
      o.SetEnableRotate(kTRUE);
   }

   // Now animate the camera
   TGLSAViewer* sav = dynamic_cast<TGLSAViewer*>(v);
   if (sav)
     sav->GetFrame()->Connect("CloseWindow()", "TTimer", &timer, "TurnOff()");
   //   timer.SetCommand("AnimateCamera()");
   timer.TurnOn();
   */

   while ( i<  10) {
      gRandom->Rannor(x, y);
      gRandom->Rannor(z,dummy);
      if ( TMath::Sqrt(x*x + y*y + z*z) < 1) {
         x = (2 * x - 1) * NucleusRadius;
         y = (2 * y - 1) * NucleusRadius;
         z = (2 * z - 1) * NucleusRadius;
         top->AddNode(proton, i, new TGeoTranslation(x, y, z));
         i++;
      }
      //      geom->CloseGeometry();
      //      geom->SetVisLevel(5);
      top->Draw("ogl");
      v->SavePicture("myAnimation.gif+");
   }
}
