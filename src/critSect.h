#include "pset.h"

class CritSect : public PSet
{
 public:
  virtual PSet* GetNewInst() {return new CritSect;}
  virtual void SanityCheck();
  virtual void PrintCurrPoss(IndStream &out, unsigned int &graphNum);
  virtual bool IsCritSect() const {return true;}
  virtual bool IsTransparent() const {return false;}
  virtual bool CanMerge(PSet *pset) const {return false;}
};
