#include "var.h"
#include <sstream>

string ModeArrayVarName(const DimVec &vec)
{
  std::stringstream name;
  name << "modes";
  DimVecConstIter iter = vec.begin();
  for(; iter != vec.end(); ++iter) {
	  name << "_" << *iter;
  }
  return name.str();
}

string IndexPairVarName(Dim dim1, Dim dim2)
{
  std::stringstream name;
  name << "indexPair_" << dim1 << "_" << dim2;
  return name.str();
}

string ModeArrayPairVarName(const DimVec &arr1, const DimVec &arr2)
{
  std::stringstream name;
  name << "modeArray";
  DimVecConstIter iter = arr1.begin();
  for(; iter != arr1.end(); ++iter) {
	  name << "_" << *iter;
  }
  name << "_";
  iter = arr2.begin();
  for(; iter != arr2.end(); ++iter) {
	  name << "_" << *iter;
  }
  return name.str();
}

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

Var::Var(Dim dim1, Dim dim2)
{
  m_type = IndexPairType;
  std::stringstream str;
  str << "c" << dim1 << "_" << dim2;
  m_compStr = str.str();
  m_pair = new std::pair<Dim,Dim>;
  m_pair->first = dim1;
  m_pair->second = dim2;
}

Var::Var(const DimVec &vec1, const DimVec &vec2)
{
  m_type = ModeArrayPairVarType;
  std::stringstream str;
  str << "d";
  DimVecConstIter iter = vec1.begin();
  for(; iter != vec1.end(); ++iter)
    str << "_" << *iter;
  iter = vec2.begin();
  for(; iter != vec2.end(); ++iter)
    str << "_" << *iter;
  m_arrPair = new std::pair<DimVec, DimVec>;
  m_arrPair->first = vec1;
  m_arrPair->second = vec2;
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
    case (IndexPairType):
      delete m_pair;
      break;
    case (ModeArrayPairVarType):
      delete m_arrPair;
      break;
    case (InvalidType) :
    default:
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
	string name = GetVarName();
	*out << "ModeArray " << name << ";\n";
	DimVecConstIter iter = m_vec->begin();
	for(; iter != m_vec->end(); ++iter) {
	  out.Indent();
	  *out << name << ".pubsh_back("
	       << *iter << ");\n";
	}
	break;
      }
    case (IndexPairType):
      {
	out.Indent();
	string name = GetVarName();
	*out << "std::pair<Index,Index> " << name << ";\n";
	break;
      }
    case (ModeArrayPairVarType):
      {
	out.Indent();
	string name = GetVarName();
	*out << "std::pair<ModeArray,ModeArray> " << name << ";\n";
	DimVecConstIter iter = m_arrPair->first.begin();
	for(; iter != m_arrPair->first.end(); ++iter) {
	  out.Indent();
	  *out << name << ".first.push_back("
	       << *iter << ");\n";
	}

	iter = m_arrPair->second.begin();
	for(; iter != m_arrPair->second.end(); ++iter) {
	  out.Indent();
	  *out << name << ".second.push_back("
	       << *iter << ");\n";
	}
      }
    case (InvalidType):
      throw;
    }
}

string Var::GetVarName() const
{
  switch (m_type) 
    {
    case (TensorVarType) :
      {
	return m_name->str();
	break;
      }
    case (ModeArrayVarType) :
      {
	return ModeArrayVarName(*(m_vec));
	break;
      }
    case (IndexPairType):
      {
	return IndexPairVarName(m_pair->first, m_pair->second);
	break;
      }
    case (ModeArrayPairVarType):
      {
	return ModeArrayPairVarName(m_arrPair->first, m_arrPair->second);
	break;
      }
    case (InvalidType):
      throw;
      return "";
    }
  throw;
}

Var::Var(const Var &var)
{
  m_type = InvalidType;
  *this = var;
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
      case (IndexPairType):
	delete m_pair;
	break;
      case (ModeArrayPairVarType):
	delete m_arrPair;
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
    case (IndexPairType):
      m_pair = new std::pair<Dim,Dim>;
      m_pair->first = rhs.m_pair->first;
      m_pair->second = rhs.m_pair->second;
      break;
    case (ModeArrayPairVarType):
      m_arrPair = new std::pair<DimVec, DimVec>;
      m_arrPair->first = rhs.m_arrPair->first;
      m_arrPair->second = rhs.m_arrPair->second;
      break;
    case (InvalidType):
      throw;
    }
  return *this;
}
