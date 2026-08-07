LoadPackage("ctbllib");;

tables := [
  "3_2.U4(3).(2^2)_{133}",
  "3_2.U4(3).2_3",
  "3_2.U4(3).2_1"
];;

Print("{\n  \"status\": \"pass\",\n  \"tables\": [\n");
for ti in [1..Length(tables)] do
  name := tables[ti];
  tbl := CharacterTable(name);
  if tbl = fail then
    Error("character table not found: ", name);
  fi;
  irr := Irr(tbl);
  orders := OrdersClassRepresentatives(tbl);
  names := ClassNames(tbl);
  pm2 := PowerMap(tbl, 2);
  if ti > 1 then Print(",\n"); fi;
  Print("    {\n");
  Print("      \"table\": \"", name, "\",\n");
  Print("      \"group_order\": \"", Size(tbl), "\",\n");
  Print("      \"class_count\": ", Length(orders), ",\n");
  Print("      \"degree_30_characters\": [\n");
  firstchar := true;
  for ci in [1..Length(irr)] do
    chi := irr[ci];
    if chi[1] = 30 then
      if not firstchar then Print(",\n"); fi;
      firstchar := false;
      Print("        {\n");
      Print("          \"character_index\": ", ci, ",\n");
      Print("          \"degree\": 30,\n");
      Print("          \"order_three_classes\": [\n");
      firstclass := true;
      has8 := false;
      for i in [1..Length(orders)] do
        if orders[i] = 3 then
          fixed := (chi[1] + chi[i] + chi[pm2[i]]) / 3;
          if not IsInt(fixed) then
            Error("fixed multiplicity is not integral");
          fi;
          if fixed = 8 then has8 := true; fi;
          if not firstclass then Print(",\n"); fi;
          firstclass := false;
          Print("            {\n");
          Print("              \"class_index\": ", i, ",\n");
          Print("              \"class_name\": \"", names[i], "\",\n");
          Print("              \"square_class_index\": ", pm2[i], ",\n");
          Print("              \"square_class_name\": \"", names[pm2[i]], "\",\n");
          Print("              \"fixed_multiplicity\": ", fixed, ",\n");
          if fixed = 8 then
            Print("              \"matches_cyclic_cubic_requirement\": true\n");
          else
            Print("              \"matches_cyclic_cubic_requirement\": false\n");
          fi;
          Print("            }");
        fi;
      od;
      Print("\n          ],\n");
      if has8 then
        Print("          \"has_fixed_rank_8_class\": true\n");
      else
        Print("          \"has_fixed_rank_8_class\": false\n");
      fi;
      Print("        }");
    fi;
  od;
  Print("\n      ]\n");
  Print("    }");
od;
Print("\n  ],\n");
Print("  \"truth_note\": \"Character-table compatibility only. A fixed multiplicity of 8 does not itself construct an integral lattice automorphism, an elliptic surface, or a rank-30 Mordell-Weil group.\"\n");
Print("}\n");
QUIT;
