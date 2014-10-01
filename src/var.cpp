/*
 This file is part of DxTer.
 DxTer is a prototype using the Design by Transformation (DxT)
 approach to program generation.
 
 Copyright (C) 2014, The University of Texas and Bryan Marker
 
 DxTer is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 DxTer is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with DxTer.  If not, see <http://www.gnu.org/licenses/>.
 */
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

string PermutationVarName(const DimVec &perm)
{
  string str = "perm";
  DimVecConstIter iter = perm.begin();
  for(; iter != perm.end(); ++iter)
    str += "_" + std::to_string(*iter);
  return str;
}

string DistEntryVecVarName(const DistEntryVec &vec)
{
  string str = "modeArrayArray";
  DistEntryVecConstIter iter = vec.begin();
  for( ; iter != vec.end(); ++iter) {
    DistEntry entry = *iter;
    DimVec vec = entry.DistEntryDims();
    DimVecConstIter iter2 = vec.begin();
    str += "__";
    for (; iter2 != vec.end(); ++iter2) {
      str += "_" + std::to_string(*iter2);
    }
  }
  return str;
}


#endif //DOTENSORS

#if DOLLDLA
string LLDLAPartVarName(const string &var, unsigned int part)
{
  std::stringstream str;
  str << var << part;
  return str.str();
}

string LLDLATransVarName(const string &var, Trans trans)
{
  switch (trans) 
    {
    case (CONJ):
      return var+"C";
    case (TRANS):
      return var+"T";
    case (CONJTRANS):
      return var+"H";
    default:
      throw;
    }
}
#endif //DOLLDLA


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
  m_compStr = "f"+m_name->str();
}
#endif

#if DOTENSORS
Var::Var(VarType type, const DimVec &vec)
{
  if (type == ModeArrayVarType) {
    m_type = ModeArrayVarType;
    m_vec = new DimVec(vec);
    std::stringstream str;
    str << "c";
    DimVecConstIter iter = m_vec->begin();
    for(; iter != m_vec->end(); ++iter)
      str << "_" << *iter;
    m_compStr = str.str();
  }
  else if (type == PermutationVarType) {
    m_type = PermutationVarType;
    m_vec = new DimVec(vec);
    std::stringstream str;
    str << "b";
    DimVecConstIter iter = m_vec->begin();
    for(; iter != m_vec->end(); ++iter)
      str << "_" << *iter;
    m_compStr = str.str();
  }
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
  m_compStr = str.str();
}


Var::Var(const DistEntryVec &vec)
{
  m_type = DistEntryVecVarType;
  m_compStr = "f " + DistEntryVecVarName(vec);
  m_entryVec = new DistEntryVec(vec);
}

#endif //DOTENSORS

#if !DOLLDLA

Var::Var(const Var &var)
{
  m_type = InvalidType;
  *this = var;
}

Var::Var(VarType type, const string &str)
{
  switch(type)
    {
    case (DirectVarDeclType):
      {
	m_type = DirectVarDeclType;
	m_varDecl = str;//new string(str);
	m_compStr = "a " + str;
	break;
      }      
#if DOTENSORS
    case (IndexArrayType):
      {
	m_type = IndexArrayType;
	m_indices = new string(str);
	//  *m_indices = indices;
	m_compStr = "e" + str;
	break;
      }
#endif //DOTENSORS
    default:
      throw;
    }
}
#endif // !DOLLDLA


#if DOLLDLA
Var::Var(const Var &var, Type dataType)
{
  m_type = InvalidType;
  *this = var;
  m_dataType = dataType;
}

Var::Var(VarType type, const string &str, Type dataType)
{
  m_dataType = dataType;
  switch(type)
    {
    case (DirectVarDeclType):
      {
	m_type = DirectVarDeclType;
        m_varDecl = str;//new string(str);
	m_compStr = "a " + str;
	break;
      }      
    default:
    {
      cout << "Error: Bad var type\n";
      throw;
    }
    }
}

Var::Var(const string &varName, unsigned int partNum, Type dataType)
{
  m_type = VarPartType;
  m_part = LLDLAPartVarName(varName, partNum);//new string(LLDLAPartVarName(varName, partNum));
  m_compStr = "g" + m_part;
  m_dataType = dataType;
}

Var::Var(const string &varName, Trans trans, Type dataType)
{
  m_type = VarTransType;
  m_transVar = new string (LLDLATransVarName(varName, trans));
  m_part = LLDLAPartVarName(varName, 0);//new string (LLDLAPartVarName(varName, 0));
  m_compStr = "f" + *m_transVar;
  m_dataType = dataType;
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
    case (PermutationVarType):
      delete m_vec;
      break;
    case (DistEntryVecVarType):
      delete m_entryVec;
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
#elif DOLLDLA
    case (VarPartType):
      //      delete m_part;
      break;
    case (VarTransType):
      delete m_transVar;
      break;
#endif
    case (DirectVarDeclType):
//      delete m_varDecl;
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
#if DOTENSORS
    case (TensorVarType) :
      {
	out.Indent();
	*out << "\t//" << m_name->PrettyStr() << endl;
	out.Indent();
	*out << "DistTensor<double> " << m_name->str() << "( "
	     << TensorDistVarName(m_name->m_type)
	     << ", g );" << endl;
	if (m_name->m_permutation.Size()) {
	  out.Indent();
	  *out << m_name->str() << ".SetLocalPerm( "
	       << PermutationVarName(m_name->m_permutation.m_permutation) 
	       << " );\n";
	}
	else {
	  Permutation defaultPerm;
	  defaultPerm.SetToDefault(m_name->m_type.m_numDims);
	  out.Indent();
	  *out << m_name->str() << ".SetLocalPerm( "
	       << PermutationVarName(defaultPerm.m_permutation) 
	       << " );\n";
	}
	break;
      }
    case (DistEntryVecVarType):
      {
	out.Indent();
	string name = GetVarName();
	*out << "std::vector<ModeArray> " << name << ";\n";
	DistEntryVecIter iter = m_entryVec->begin();
	for(; iter != m_entryVec->end(); ++iter) {
	  DistEntry entry = *iter;
	  out.Indent();
	  *out << name << ".push_back(" <<
	    ModeArrayVarName(entry.DistEntryDims()) 
	       << ");\n";
	}
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
    case (PermutationVarType):      
      {
	out.Indent();
	string name = GetVarName();
	*out << "Permutation " << name << ";\n";
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
	*out << "std::pair<Index,Index> " << name 
	     << "(" << m_pair->first << "," << m_pair->second << ");\n";
	
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
	*out <<"]";
	if (!m_distType->m_notReped.IsStar()) {
	  *out << "|(";
	  DimVec vec = m_distType->m_notReped.DistEntryDims();
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
	*out << "\");\n";
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
#elif DOLLDLA
    case (VarPartType):
      {
	out.Indent();
	if (m_dataType == REAL_SINGLE) {
	  *out << "float *" << m_part << ";\n";
	} else if (m_dataType == REAL_DOUBLE) {
	  *out << "double *" << m_part << ";\n";
	} else {
	  cout << "ERROR: Var " << m_part << " has invalid m_dataType\n";
	  throw;
	}
	break;
      }
    case (VarTransType):
      {
	out.Indent();
	if (m_dataType == REAL_SINGLE) {
	  *out << "float *" << m_part << ";\n";
	} else if (m_dataType == REAL_DOUBLE) {
	  *out << "double *" << m_part << ";\n";
	} else {
	  cout << "ERROR: Var " << m_part << " has invalid m_dataType\n";
	  throw;
	}
	break;
      }
#endif
    case (DirectVarDeclType):
      {
	out.Indent();
//	cout << "DirectVarDeclType " << *m_varDecl << endl;
//	m_varDecl = new string("Nope\n");
	*out << m_varDecl << endl;
	break;
      }
    case (InvalidType):
      {
	cout << "Error: Invalid var type\n";
	throw;
      }
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
    case (DistEntryVecVarType) :
      {
	return DistEntryVecVarName(*m_entryVec);
	break;
      }
    case (ModeArrayVarType) :
      {
	return ModeArrayVarName(*(m_vec));
	break;
      }
    case (PermutationVarType) :
      {
	return PermutationVarName(*(m_vec));
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
#elif DOLLDLA
    case (VarPartType):
      {
	return m_part;
	break;
      }
    case (VarTransType):
      {
	return *m_transVar;
	break;
      }
#endif
    case (DirectVarDeclType):
      throw;
    case (InvalidType):
      throw;
      return "";
    }
  throw;
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
      case (PermutationVarType) :
	delete m_vec;
	break;
      case (DistEntryVecVarType):
	delete m_entryVec;
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
#elif DOLLDLA
      case (VarPartType):
	//	delete m_part;
	break;
      case (VarTransType):
	delete m_transVar;
	break;
#endif
      case (DirectVarDeclType):
	//	delete m_varDecl;
	break;
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
    case (PermutationVarType) :
      m_vec = new DimVec (*(rhs.m_vec));
      break;
    case (DistEntryVecVarType):
      m_entryVec = new DistEntryVec(*(rhs.m_entryVec));
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
#elif DOLLDLA
    case (VarPartType):
      m_part = rhs.m_part;//new string(*(rhs.m_part));
      break;
    case (VarTransType):
      m_transVar = new string (*(rhs.m_transVar));
      break;
#endif
    case (DirectVarDeclType):
      m_varDecl = rhs.m_varDecl;//new string ((rhs.m_varDecl));
      break;
    case (InvalidType):
      throw;
    }
  return *this;
}
