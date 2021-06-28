#ifndef MODELCONTROLLERINTERFACE_H
#define MODELCONTROLLERINTERFACE_H


class ModelControllerInterface
{
public:
    virtual void Prepare(void) = 0;
    virtual bool Step(void) = 0;
    virtual void RequestAbort(void) = 0;
};


#endif // MODELCONTROLLERINTERFACE_H
