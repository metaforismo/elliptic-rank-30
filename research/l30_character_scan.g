LoadPackage("ctbllib");;

tables := [
  "3_2.U4(3).(2^2)_{133}",
  "3_2.U4(3).2_3",
  "3_2.U4(3).2_1"
];;

# Emit a compact TSV stream containing integers only.  We record both
# absolutely irreducible complex characters and Galois-rationalized character
# orbits.  The 30-dimensional lattice representation may be rationally
# irreducible while splitting over a cyclotomic field.
for ti in [1..Length(tables)] do
  tbl := CharacterTable(tables[ti]);
  if tbl = fail then
    Error("character table not found: ", tables[ti]);
  fi;
  irr := Irr(tbl);
  irrvals := List(irr, ValuesOfClassFunction);
  ratchars := RationalizedMat(irrvals);
  orders := OrdersClassRepresentatives(tbl);
  pm2 := PowerMap(tbl, 2);
  Print("TABLE|", ti, "|", Size(tbl), "|", Length(orders), "\n");
  for ci in [1..Length(irr)] do
    chi := irr[ci];
    if chi[1] = 30 then
      Print("CHAR|", ti, "|", ci, "\n");
      for i in [1..Length(orders)] do
        if orders[i] = 3 then
          fixed := (chi[1] + chi[i] + chi[pm2[i]]) / 3;
          if not IsInt(fixed) then Error("complex fixed multiplicity nonintegral"); fi;
          Print("CLASS|", ti, "|", ci, "|", i, "|", pm2[i], "|", fixed, "\n");
        fi;
      od;
    fi;
  od;
  for ri in [1..Length(ratchars)] do
    chi := ratchars[ri];
    if chi[1] <= 60 then
      Print("RAT|", ti, "|", ri, "|", chi[1], "\n");
      if chi[1] = 30 then
        for i in [1..Length(orders)] do
          if orders[i] = 3 then
            fixed := (chi[1] + chi[i] + chi[pm2[i]]) / 3;
            if not IsInt(fixed) then Error("rational fixed multiplicity nonintegral"); fi;
            Print("RCLASS|", ti, "|", ri, "|", i, "|", pm2[i], "|", fixed, "\n");
          fi;
        od;
      fi;
    fi;
  od;
od;
QUIT;
