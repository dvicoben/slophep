#ifndef SL_Decay_Config_h
#define SL_Decay_Config_h

#include<string>

class BR_mapping
{
private:
    std::string dacay_tag;
    double BR_value;
    
public:
    std::string get_BR_tag(void) const;
    double get_BR(void) const;
    BR_mapping(std::string s, double val);
    ~BR_mapping();
};

#endif // SL_Decay_Config_h
