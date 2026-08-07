LoadPackage("atlasrep");;
SetUserPreference("AtlasRep", "RemoteAccess", true);;

names := [
  "6_1.U4(3).2_2",
  "3_1.U4(3).2_2",
  "12_1.U4(3).2_2"
];;

for ni in [1..Length(names)] do
  name := names[ni];
  infos := AllAtlasGeneratingSetInfos(name);
  Print("GROUP|", ni, "|", name, "|", Length(infos), "\n");
  for i in [1..Length(infos)] do
    info := infos[i];
    dim := -1;
    p := -1;
    ring := "";
    charname := "";
    repname := "";
    if IsBound(info.dim) then dim := info.dim; fi;
    if IsBound(info.p) then p := info.p; fi;
    if IsBound(info.ring) then ring := String(info.ring); fi;
    if IsBound(info.charactername) then charname := info.charactername; fi;
    if IsBound(info.repname) then repname := info.repname; fi;
    if dim <= 60 or dim = -1 then
      Print("INFO|", ni, "|", i, "|", dim, "|", p, "|", charname, "|", repname, "|", ring, "\n");
    fi;
  od;
od;
QUIT;
