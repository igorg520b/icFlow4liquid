#ifndef MODELCONTROLLERINTERFACE_H
#define MODELCONTROLLERINTERFACE_H


class ModelControllerInterface
{
public:
    virtual void Prepare() = 0;
    virtual bool Step() = 0;
    virtual void RequestAbort() = 0;
};


#endif // MODELCONTROLLERINTERFACE_H
