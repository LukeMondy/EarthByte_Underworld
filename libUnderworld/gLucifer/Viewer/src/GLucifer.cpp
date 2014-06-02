/********************************************************************************************************************** 
 * Attempt to get gLucifer running in OmegaLib
 *********************************************************************************************************************/
#ifdef USE_OMEGALIB

#include <omega.h>
#include <omegaGl.h>
#include "ViewerApp.h"
#include "GLuciferViewer.h"

std::vector<std::string> args;

using namespace omega;

class GLuciferApplication;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class GLuciferRenderPass: public RenderPass, ViewerApp
{
public:
	GLuciferRenderPass(Renderer* client, GLuciferApplication* app, OpenGLViewer* viewer): RenderPass(client, "GLuciferRenderPass"), app(app), ViewerApp(viewer) {}
	virtual void initialize();
	virtual void render(Renderer* client, const DrawContext& context);

   // Virtual functions for interactivity (from ViewerApp/ApplicationInterface)
   virtual bool mouseMove(int x, int y) {}
   virtual bool mousePress(MouseButton btn, bool down, int x, int y) {}
   virtual bool mouseScroll(int scroll) {}
   virtual bool keyPress(unsigned char key, int x, int y) {}

private:
	GLuciferApplication* app;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class GLuciferApplication: public EngineModule
{
public:
  OpenGLViewer* viewer;
  GLuciferViewer* glapp;
  bool redisplay = true;
  int argc;
  char** argv;

	GLuciferApplication(): EngineModule("GLuciferApplication") { enableSharedData(); }

	virtual void initializeRenderer(Renderer* r) 
	{ 
    viewer = new OpenGLViewer(false, false);
		r->addRenderPass(new GLuciferRenderPass(r, this, viewer));
	}

	float setArgs(int argc, char** argv) {this->argc = argc; this->argv = argv;}

	float getYaw() { return myYaw; }
	float getPitch() { return myPitch; }

	virtual void handleEvent(const Event& evt);
	virtual void commitSharedData(SharedOStream& out);
	virtual void updateSharedData(SharedIStream& in);

private:
	float myYaw;
	float myPitch;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GLuciferRenderPass::initialize()
{
	RenderPass::initialize();

  //Init fractal app
		DisplaySystem* ds = app->getEngine()->getDisplaySystem();
		Vector2i resolution = ds->getCanvasSize();
  //Fake the arg list from vector of args
  int argc = args.size()+1;
  char* argv[argc];
  argv[0] = (char*)malloc(20);
  strcpy(argv[0], "./gLucifer");
  for (int i=1; i<argc; i++)
    argv[i] = (char*)args[i-1].c_str();
  //Create the app
  app->glapp = new GLuciferViewer(args, viewer, resolution[0], resolution[1]);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GLuciferRenderPass::render(Renderer* client, const DrawContext& context)
{
	if(context.task == DrawContext::SceneDrawTask)
	{
		client->getRenderer()->beginDraw3D(context);

    if (!viewer->isopen)
    {
      //Load vis data for first window
      if (!app->glapp->loadWindow(0, -1, true))
         abort_program("Model file load error, no window data\n");

      DisplaySystem* ds = app->getEngine()->getDisplaySystem();
      Vector2i resolution = ds->getCanvasSize();
      viewer->open(resolution[0], resolution[1]);
    }
    //if (app->redisplay)
    {
      viewer->display();
      app->redisplay = false;
      //printf("display\n");
    }

		client->getRenderer()->endDraw();
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GLuciferApplication::handleEvent(const Event& evt)
{
  //printf(". %d %d %d\n", evt.getType(), evt.getServiceType(), evt.getFlags());
	if(evt.getServiceType() == Service::Pointer)
	{
    int x = evt.getPosition().x();
    int y = evt.getPosition().y();
    int flags = evt.getFlags();
    MouseButton button = NoButton;
    if (flags & 1)
      button = LeftButton;
    else if (flags & 2)
      button = RightButton;
    else if (flags & 4)
      button = MiddleButton;

    switch (evt.getType())
    {
      case Event::Down:
        //printf("%d %d\n", button, flags);
        if (button <= RightButton) viewer->mouseState ^= (int)pow(2, (int)button);
        viewer->mousePress(button, true, x, y);
          redisplay = true;
          break;
      case Event::Up:
        viewer->mouseState = 0;
        viewer->mousePress(button, false, x, y);
          break;
      case Event::Zoom:
        viewer->mouseScroll(evt.getExtraDataInt(0));
         break;
      case Event::Move:
         if (viewer->mouseState)
         {
            viewer->mouseMove(x, y);
            //redisplay = true;
         }
          break;
      default:
        printf("? %d\n", evt.getType());
    }



		/*/ Normalize the mouse position using the total display resolution, 
		// then multiply to get 180 degree rotations
		DisplaySystem* ds = getEngine()->getDisplaySystem();
		Vector2i resolution = ds->getCanvasSize();
		myYaw = (evt.getPosition().x() / resolution[0]) * 180;
		myPitch = (evt.getPosition().y() / resolution[1]) * 180;
    */
	}
  else if(evt.getServiceType() == Service::Keyboard)
  {
    int x = evt.getPosition().x();
    int y = evt.getPosition().y();
    int key = evt.getSourceId();
    if (evt.isKeyDown(key))
    {
      if (key > 255)
      {
      //printf("Key %d %d\n", key, evt.getFlags());
        if (key == 262) key = KEY_UP;
        else if (key == 264) key = KEY_DOWN;
        else if (key == 261) key = KEY_LEFT;
        else if (key == 263) key = KEY_RIGHT;
        else if (key == 265) key = KEY_PAGEUP;
        else if (key == 266) key = KEY_PAGEDOWN;
        else if (key == 260) key = KEY_HOME;
        else if (key == 267) key = KEY_END;
      }
      viewer->keyPress(key, x, y);
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GLuciferApplication::commitSharedData(SharedOStream& out)
{
	out << myYaw << myPitch;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GLuciferApplication::updateSharedData(SharedIStream& in)
{
 	in >> myYaw >> myPitch;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ApplicationBase entry point
int main(int argc, char** argv)
{
	Application<GLuciferApplication> app("gLucifer");
  oargs().setStringVector("gLucifer", "gLucifer Arguments", args);
  return omain(app, argc, argv);
}

#endif
