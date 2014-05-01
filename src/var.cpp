#include "var.h"
#include <sstream>

Var::Var(const Name &name)
{
  m_type = TensorVarType;
  m_name = new Name(name);
  m_compStr = "a"+m_name->str();
}


Var::Var(const DimVec &vec)
{
  m_type = ModeArrayVarType;
  m_vec = new DimVec(vec);
  std::stringstream str;
  str << "b";
  DimVecConstIter iter = m_vec->begin();
  for(; iter != m_vec->end(); ++iter)
    str << "_" << *iter;
  m_compStr = str.str();
}

Var::~Var()
{
  switch (m_type) 
    {
    case (TensorVarType) :
      delete m_name;
      break;
    case (ModeArrayVarType) :
      delete m_vec;
      break;
    case (InvalidType) :
      throw;
    }
}


void Var::PrintDecl(IndStream &out) const
{
  switch (m_type) 
    {
    case (TensorVarType) :
      {
	out.Indent();
	*out << "\t//" << m_name->PrettyStr() << endl;
	out.Indent();
	*out << "DistTensor<double> " << m_name->str() << "(shape, dist, indices, g);" << endl;
	break;
      }
    case (ModeArrayVarType) :
      {
	out.Indent();
	std::stringstream name;
	name << "modes";
	DimVecConstIter iter = m_vec->begin();
	for(; iter != m_vec->end(); ++iter) {
	  name << "_" << *iter;
	}
	*out << "ModeArray " << name.str() << ";\n";
	out.Indent();
	iter = m_vec->begin();
	for(; iter != m_vec->end(); ++iter) {
	  *out << name.str() << ".push_back("
	       << *iter << ");\n";
	}
	break;
      }
    case (InvalidType):
      throw;
    }
}


Var::Var(const Var &var)
{
  m_type = var.m_type;
  m_compStr = var.m_compStr;
  switch (m_type) 
    {
    case (TensorVarType):
      m_name = new Name(*(var.m_name));
      break;
    case (ModeArrayVarType):
      m_vec = new DimVec (*(var.m_vec));
      break;
    case (InvalidType):
      throw;
    }
}

Var& Var::operator=(const Var &rhs)
{
  if (m_type != InvalidType) {
    switch (m_type) 
      {
      case (TensorVarType) :
	delete m_name;
	break;
      case (ModeArrayVarType) :
	delete m_vec;
	break;
      case (InvalidType):
	throw;
      }    
  }
  m_type = rhs.m_type;
  m_compStr = rhs.m_compStr;
  switch (m_type) 
    {
    case (TensorVarType):
      m_name = new Name(*(rhs.m_name));
      break;
    case (ModeArrayVarType):
      m_vec = new DimVec (*(rhs.m_vec));
      break;
    case (InvalidType):
      throw;
    }
  return *this;
}
