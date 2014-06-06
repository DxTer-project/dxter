#include "var.h"
#include <sstream>

#if DOTENSORS
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

string TensorDistVarName(const DistType &type)
{
  return "dist_" + type.str();
}

string IndexArrayVarName(const string &indices)
{
  return "indices_" + indices;
}
#endif

#if DOTENSORS
Var::Var(const DistType &type)
{
  m_type = TensorDistVarType;
  m_distType = new DistType(type);
  m_compStr = "a"+type.QuickStr();
}
#endif

#if DOTENSORS
Var::Var(const Name &name)
{
  m_type = TensorVarType;
  m_name = new Name(name);
  m_compStr = "b"+m_name->str();
}
#endif

#if DOTENSORS
Var::Var(const DimVec &vec)
{
  m_type = ModeArrayVarType;
  m_vec = new DimVec(vec);
  std::stringstream str;
  str << "c";
  DimVecConstIter iter = m_vec->begin();
  for(; iter != m_vec->end(); ++iter)
    str << "_" << *iter;
  m_compStr = str.str();
}


Var::Var(Dim dim1, Dim dim2)
{
  m_type = IndexPairType;
  std::stringstream str;
  str << "d" << dim1 << "_" << dim2;
  m_compStr = str.str();
  m_pair = new std::pair<Dim,Dim>;
  m_pair->first = dim1;
  m_pair->second = dim2;
}


Var::Var(const DimVec &vec1, const DimVec &vec2)
{
  m_type = ModeArrayPairVarType;
  std::stringstream str;
  str << "e";
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

Var::Var(const string &indices)
{
  m_type = IndexArrayType;
  m_indices = new string(indices);
  //  *m_indices = indices;
  m_compStr = "e" + indices;
}

#endif

Var::~Var()
{
  switch (m_type) 
    {
#if DOTENSORS
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
    case (TensorDistVarType):
      delete m_distType;
      break;
    case (IndexArrayType):
      delete m_indices;
      break;
#endif
    case (InvalidType) :
    default:
      throw;
    }
}


void Var::PrintDecl(IndStream &out) const
{
  switch (m_type) 
    {
#if DOTENSORS
    case (TensorVarType) :
      {
	out.Indent();
	*out << "\t//" << m_name->PrettyStr() << endl;
	out.Indent();
	*out << "DistTensor<double> " << m_name->str() << "( "
	     << TensorDistVarName(m_name->m_type)
	     << ", g );" << endl;
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
	  *out << name << ".push_back("
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
	break;
      }
    case (TensorDistVarType):
      {
	out.Indent();
	*out << "TensorDistribution " << GetVarName() << " = "
	     << "tmen::StringToTensorDist(\"[";
	for(Dim dim = 0; dim < m_distType->m_numDims; ++dim) {
	  DimVec vec = m_distType->m_dists[dim].DistEntryDims();
	  if (dim)
	    *out << ",";
	  *out << "(";
	  bool start = true;
	  DimVecIter iter = vec.begin();
	  for( ; iter != vec.end(); ++iter) {
	    if (!start) {
	      *out << ",";
	    }
	    else
	      start = false;
	    *out << *iter;
	  }
	  *out << ")";
	}
	*out <<"]\");\n";
	break;
      }
    case (IndexArrayType):
      {
	out.Indent();
	string name = GetVarName();
	*out << "IndexArray " << name << "( " << m_indices->size() << " );\n";
	string::iterator iter = m_indices->begin();
	for(unsigned int i = 0; iter != m_indices->end(); ++iter,++i) {
	  out.Indent();
	  *out << name << "[" << i <<"] = '" << *iter << "';\n";
	}
	break;
      }
#endif
    case (InvalidType):
      throw;
    }
}

string Var::GetVarName() const
{
  switch (m_type) 
    {
#if DOTENSORS
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
    case (TensorDistVarType):
      {
	return TensorDistVarName(*m_distType);
	break;
      }
    case (IndexArrayType):
      {
	return IndexArrayVarName(*m_indices);
	break;
      }
#endif
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
#if DOTENSORS
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
      case (TensorDistVarType):
	delete m_distType;
	break;
      case (IndexArrayType):
	delete m_indices;
	break;
#endif
      case (InvalidType):
	throw;
      }    
  }
  m_type = rhs.m_type;
  m_compStr = rhs.m_compStr;
  switch (m_type) 
    {
#if DOTENSORS
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
    case (TensorDistVarType):
      m_distType = new DistType(*rhs.m_distType);
      break;
    case (IndexArrayType):
      m_indices = new std::string(*(rhs.m_indices));
      break;
#endif
    case (InvalidType):
      throw;
    }
  return *this;
}
