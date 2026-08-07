LoadPackage("ctbllib");;
SizeScreen([10000, 10000]);;

names := AllCharacterTableNames(
  Identifier,
  x -> PositionSublist(x, "U4(3)") <> fail
);;
relevantOrders := [19595520, 39191040, 78382080];;

for ni in [1..Length(names)] do
  name := names[ni];
  tbl := CharacterTable(name);
  if tbl <> fail and Size(tbl) in relevantOrders then
    irr := Irr(tbl);
    irrvals := List(irr, ValuesOfClassFunction);
    rat := RationalizedMat(irrvals);
    orders := OrdersClassRepresentatives(tbl);
    pm2 := PowerMap(tbl, 2);
    indicators := Indicator(tbl, irr, 2);
    Print("TABLE|", ni, "|", name, "|", Size(tbl), "|", Length(orders), "\n");

    for ci in [1..Length(irr)] do
      chi := irr[ci];
      if chi[1] = 15 or chi[1] = 30 then
        Print("COMPLEX|", ni, "|", ci, "|", chi[1], "|", indicators[ci], "\n");
        for i in [1..Length(orders)] do
          if orders[i] = 3 then
            fixed := (chi[1] + chi[i] + chi[pm2[i]]) / 3;
            if not IsInt(fixed) then Error("nonintegral fixed multiplicity"); fi;
            Print("CCLASS|", ni, "|", ci, "|", i, "|", pm2[i], "|", fixed, "\n");
          fi;
        od;
      fi;
    od;

    for ri in [1..Length(rat)] do
      chi := rat[ri];
      if chi[1] = 15 or chi[1] = 30 then
        Print("RATIONAL|", ni, "|", ri, "|", chi[1], "\n");
        for i in [1..Length(orders)] do
          if orders[i] = 3 then
            fixed := (chi[1] + chi[i] + chi[pm2[i]]) / 3;
            if not IsInt(fixed) then Error("nonintegral rationalized fixed multiplicity"); fi;
            Print("RCLASS|", ni, "|", ri, "|", i, "|", pm2[i], "|", fixed, "\n");
          fi;
        od;
      fi;
    od;
  fi;
od;
QUIT;
