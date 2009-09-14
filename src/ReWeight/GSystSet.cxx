//____________________________________________________________________________
/*
 Copyright (c) 2003-2009, GENIE Neutrino MC Generator Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Costas Andreopoulos <costas.andreopoulos \at stfc.ac.uk>
         STFC, Rutherford Appleton Laboratory 

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :
 @ Aug 01, 2009 - CA
   Was adapted from Jim's and Costas' T2K-specific GENIE reweighting code.
   First included in v2.5.1.

*/
//____________________________________________________________________________

#include <iostream>

#include "ReWeight/GSystSet.h"

using std::cerr;
using std::cout;
using std::endl;

using namespace genie;
using namespace genie::rew;

//_______________________________________________________________________________________
GSystSet::GSystSet()
{

}
//_______________________________________________________________________________________
GSystSet::GSystSet(const GSystSet & syst)
{
  this->Copy(syst);
}
//_______________________________________________________________________________________
GSystSet::~GSystSet()
{
  fSystematics.clear();
}
//_______________________________________________________________________________________
void GSystSet::Include(GSyst_t syst)
{
  if(syst == kSystNull) return;

  GSystInfo * syst_info = new GSystInfo;
  fSystematics.insert( map<GSyst_t, GSystInfo*>::value_type(syst, syst_info) );     
}
//_______________________________________________________________________________________
int GSystSet::NIncluded(void) const
{
  return fSystematics.size();
}
//_______________________________________________________________________________________
bool GSystSet::IsIncluded(GSyst_t syst) const
{
  return (fSystematics.find(syst) != fSystematics.end());
}
//_______________________________________________________________________________________
vector<genie::rew::GSyst_t> GSystSet::AllIncluded(void)
{
  vector<GSyst_t> svec;

  map<GSyst_t, GSystInfo*>::const_iterator it = fSystematics.begin();
  for( ; it != fSystematics.end(); ++it) {
    GSyst_t syst = it->first;
    svec.push_back(syst);
  }
  return svec;
}
//_______________________________________________________________________________________
void GSystSet::Remove(GSyst_t syst)
{
  fSystematics.erase(syst);
}
//_______________________________________________________________________________________
double GSystSet::CurValue(GSyst_t syst) const
{
  if ( this->IsIncluded(syst) ) return fSystematics.find(syst)->second->CurValue;
  else return 0.;
}
//_______________________________________________________________________________________
double GSystSet::DefValue(GSyst_t syst) const
{
  if ( this->IsIncluded(syst) ) return fSystematics.find(syst)->second->DefValue;
  else return 0.;
}
//_______________________________________________________________________________________
double GSystSet::MinValue(GSyst_t syst) const
{
  if ( this->IsIncluded(syst) ) return fSystematics.find(syst)->second->MinValue;
  else return 0.;
}
//_______________________________________________________________________________________
double GSystSet::MaxValue(GSyst_t syst) const
{
  if ( this->IsIncluded(syst) ) return fSystematics.find(syst)->second->MaxValue;
  else return 0.;
}
//_______________________________________________________________________________________
void GSystSet::SetCurValue(GSyst_t syst, double val)
{
  if ( this->IsIncluded(syst) ) {
    fSystematics[syst]->CurValue = val;
  }
}
//_______________________________________________________________________________________
void GSystSet::SetDefValue(GSyst_t syst, double val)
{
  if ( this->IsIncluded(syst) ) {
    fSystematics[syst]->DefValue = val;
  }
}
//_______________________________________________________________________________________
void GSystSet::SetRange(GSyst_t syst, double min, double max)
{
  if ( this->IsIncluded(syst) ) {
    fSystematics[syst]->MinValue = min;
    fSystematics[syst]->MaxValue = max;
  }
}
//_______________________________________________________________________________________
void GSystSet::PrintSummary(void)
{
  cout << "Have included " << this->NIncluded() << " systematics:" << endl;  
				    
  vector<genie::rew::GSyst_t> svec = this->AllIncluded();
  GSyst syst_obj;

  vector<genie::rew::GSyst_t>::const_iterator it = svec.begin();
  for( ; it != svec.end(); ++it) {
    GSyst_t syst = *it;
    
    cout << " |--> " << syst_obj.AsString(syst) << endl;
  }

}
//_______________________________________________________________________________________
void GSystSet::Copy(const GSystSet & syst_set)
{
  return fSystematics.clear();

  map<GSyst_t, GSystInfo*>::const_iterator it = syst_set.fSystematics.begin();
  for( ; it != syst_set.fSystematics.end(); ++it) {
    GSyst_t     syst       = it->first;
    GSystInfo * syst_info  = it->second;

    double cur = syst_info->CurValue;
    double def = syst_info->DefValue;
    double min = syst_info->MinValue;
    double max = syst_info->MaxValue;

    this->Include     (syst);
    this->SetCurValue (syst, cur);
    this->SetDefValue (syst, def);
    this->SetRange    (syst, min, max);
  }
}
//_______________________________________________________________________________________