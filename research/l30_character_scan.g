LoadPackage("ctbllib");;

tables := [
  "3_2.U4(3).(2^2)_{133}",
  "3_2.U4(3).2_3",
  "3_2.U4(3).2_1"
];;

# Emit a compact, machine-readable TSV stream containing integers only.
# The GitHub workflow adds the stable table names and serializes JSON in Python.
for ti in [1..Length(tables)] do
  tbl := CharacterTable(tables[ti]);
  if tbl = fail then
    Error("character table not found: ", tables[ti]);
  fi;
  irr := Irr(tbl);
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
          if not IsInt(fixed) then
            Error("fixed multiplicity is not integral");
          fi;
          Print("CLASS|", ti, "|", ci, "|", i, "|", pm2[i], "|", fixed, "\n");
        fi;
      od;
    fi;
  od;
od;
QUIT;
