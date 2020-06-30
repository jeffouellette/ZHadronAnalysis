#ifndef __Trigger_h__
#define __Trigger_h__

/**
 * Implements a trigger class inspired by a linked list design. A trigger points at its reference trigger, which is used to define its efficiency factor.
 * Author: Jeff Ouellette
 * Dated: 4/25/2018
 */

#include <string>

using namespace std;

/**
 * Trigger stores information about a trigger, including a jet momenta range, pseudorapidity interval,
 * name, and a (unique) branching index for event analysis.
 */
struct Trigger {
    
  public:
    string name = "";

    int minPt = 0;
    int maxPt = 0;
    int thresPt = 0;
    double lowerEta = 0;
    double upperEta = 0;
    int lowerRunNumber = 0;
    int upperRunNumber = 0;
    int index = 0;
    bool disabled = false;
    bool isBootstrapped = false;
    Trigger* referenceTrigger = NULL;

    bool trigBool = false;
    float trigPrescale = -1;

    Trigger (const string _name, const int _thresPt, const double _lowerEta, const double _upperEta, const int _lRN=0, const int _uRN=10000000, const int _minPt=0);
    Trigger (const string _name, const int _thresPt, const double _lowerEta, const double _upperEta, const bool _disabled, const int _lRN=0, const int _uRN=10000000, const int _minPt=0);
    Trigger (const Trigger* t);

    ~Trigger () {}

};

#endif
