(* Fichier ttgraphique.ml : librairie graphique adaptée à la représentation de graphes sous forme de tableaux, listes, etc.*)
#open "graphics";;
open_graph "";;

let trace_cadre tx ty pasx pasy gridbool = (*trace un cadre gradué, les tailles sont en px*)
  set_color 0x000000;
  let sx =size_x() and sy = size_y() in
  moveto 20 20;
  lineto 20 (20+ty);
  moveto 20 20;
  lineto (20+tx) 20;
  for i = 1 to ty/pasy do
    moveto 20 (20+i*pasy);
    lineto 24 (20+i*pasy);
    if gridbool then
    begin
      set_color 0xbdbdbd;
      lineto (20+tx) (20+i*pasy);
      set_color 0x000000;
    end;
  done;
  for i = 1 to tx/pasx do
    moveto (20+i*pasx) 20;
    lineto (20+i*pasx) 24;
    if gridbool then
    begin
      set_color 0xbdbdbd;
      lineto (20+i*pasx) (20+ty);
      set_color 0x000000;
    end;
  done;
;;

let conversion data pas =
  let rec c_aux dat posi= (*on transforme un ensemble de données (ordonnées ? à voir) en un ensemble de couples de points ordonné selon une abscisse...*)
    match dat with
      |[]->[];
      |a::q->(posi,a)::(c_aux q (posi+pas));
  in
  c_aux data 0;
;;

let rec trace_dot points mult = match points with
  |[]->();
  |[(x,y)]->();
  |(x1,y1)::(x2,y2)::q->
    let (y1t,y2t) = (int_of_float (mult*.y1),int_of_float (mult*.y2)) in
    moveto (x1+20) (y1t+20);
    lineto (x2+20) (y2t+20);
    trace_dot ((x2,y2)::q) mult;;

(* Routine test graphique*)
let (N,P) = initpmc ();;
let bas = "mettre chose ici";;(* apprentissage, qui a des effets de bord et renvoie la liste des erreurs globales *)

let dots = conversion bas 1;;

clear_graph();;

trace_cadre 1000 600 2 30 false;;

set_color 0x000000;;

trace_dot dots 600.;;




