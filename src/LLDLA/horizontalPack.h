/*
    This file is part of DxTer.
    DxTer is a prototype using the Design by Transformation (DxT)
    approach to program generation.

    Copyright (C) 2015, The University of Texas and Bryan Marker

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

#if DOLLDLA

class HorizontalPack : public Pack {
 public:
  explicit Pack(Layer layer);
  static Node* BlankInst() { return new HorizontalPack(ABSLAYER); }
  virtual Node* GetNewInst() { return BlankInst(); }

  virtual NodeType GetType() const { return "HorizontalPack"; }
  virtual ClassType GetNodeClass() const { return GetClass(); }
  static ClassType GetClass() { return "HorizontalPack"; }

  virtual void Prop();
};
#endif // DOLLDLA
