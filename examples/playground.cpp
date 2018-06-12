
#include <iostream>
#include <pfsoft>

using namespace pfsoft;

template <template <class> class CreationPolicy>
class WidgetManager : public CreationPolicy<double>
{
    void DoSomething()
    {
        
    }
};

int main()
{
    PFSOFT_STATIC_CHECK(sizeof(int) < sizeof(double), Destination_Type_Too_Narrow);
    WidgetManager<OpNewCreator> MyWidgetManager;
    
    Chunk c;
}
