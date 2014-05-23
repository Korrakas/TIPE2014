(* Ce fichier contient la totalité des fonctions de base, modifiées pour pouvoir faire des expériences rapidement
sur la base du PMC sans inertie*)
(*TODO : ajout d'un biais et d'une regle de gestion du biais*)

random__init (int_of_float (sys__time()));;

let config = [|3;3;7;2|];; (* dans config, on pense à ajouter le biais dans chaque couche : le biais est en j=0 de manière SYSTEMATIQUE !*)

let tailleH = vect_length config;;

let print_floatvect v =
	print_string "[|";
	for i = 0 to vect_length v - 1 do
		print_float v.(i);
		print_string "; ";
	done;
	print_string "|]";
;;

let initpmc () =
	let neur, poids = make_vect tailleH [||], make_vect tailleH [||]
	in
	let fillp pij t =
		for k = 1 to t-1 do (*pas de lien vers le biais*)
			pij.(k) <- random__float 0.1 -. 0.05; (* on remplit le tableau des poids synaptiques aléatoirement entre 0.5 et 1.*)
		done;
	in
	for i = 0 to tailleH - 1 do
		neur.(i) <- make_vect config.(i) 0.;
		poids.(i) <- make_vect config.(i) [||];
		neur.(i).(0) <- -1.;
		if i <> tailleH - 1 then (*pas de liens externes pour la couche de sortie*) 
		for j = 0 to config.(i) - 1 do
			poids.(i).(j) <- make_vect config.(i+1) 0.;
			fillp poids.(i).(j) config.(i+1);
		done;
	done;
	neur, poids;
;;



let act x = 1./.(1. +. exp(-.x));;

let actprime x = exp(x)/.(power (exp(x) +. 1.) 2.);;

let actinv x = 0. -. log(1./.x -. 1.);;

exception taille_incompatible;;

let entree neur tab =
	let tailletab = vect_length tab in
	if tailletab <> config.(0) then raise taille_incompatible;
	for j = 0 to config.(0)-1 do
		neur.(0).(j) <- tab.(j);
	done;
;;

let propagation (neur, poids) =
let temp = ref 0. in
	for i = 1 to tailleH -1 do
		for j = 1 to config.(i) - 1 do (*on commence les modifs à j=1, car j=0 est réservé au biais.*)
			temp := 0.;
			for h = 0 to config.(i-1) -1 do (*on commence à h=0, le biais entre bien dans le calcul toutefois.*)
				temp := !temp +. poids.(i-1).(h).(j)*.neur.(i-1).(h) ;
			done;
			neur.(i).(j) <- act(!temp);
		done;
	done;
;;

let erreurglobale neur sortie =
	let taillesortie = vect_length sortie and er =  ref 0. in
	if taillesortie<>config.(tailleH-1) then raise taille_incompatible;
	for j = 0 to config.(tailleH-1)-1 do
		er := !er +. power (sortie.(j) -. neur.(tailleH-1).(j)) 2.;
	done;
	!er;
;;

let creergradient (neur,poids) sortie =
	let grad = map_vect copy_vect neur and temp = ref 0. in
	for j = 0 to config.(tailleH-1)-1 do
		grad.(tailleH - 1).(j) <- sortie.(j) -. neur.(tailleH -1).(j)
	done;
	for i = tailleH -2 downto 1 do
		for j = 0 to config.(i)-1 do
			temp := 0.;
			for h = 1 to config.(i+1)-1 do (*rétropropagation : on oublie le gradient de (0), le biais, non ?*)
				temp:= !temp +. grad.(i+1).(h)*.poids.(i).(j).(h)
			done;
			grad.(i).(j) <- !temp;
		done;
	done;
	grad;
;;


let modifpoids (neur,poids) grad erg nu=
	for i = 0 to tailleH - 2 do
		for j = 0 to config.(i) - 1 do
			for h = 1 to config.(i+1) -1 do (*pas besoin de toucher le "poids vers le biais" qui ne doit de toute façon pas exister (0 ?) *)
				poids.(i).(j).(h) <- poids.(i).(j).(h) +. nu*.grad.(i+1).(h)*.(actprime erg)*.neur.(i).(j); 
			done;
		done;
	done;
;;

let resolution = 5;; (*dans la liste ne sont intégrées que les données tous les 'resolution' pas*)

let apprentissage1 (N,P) (in_data,out_data) iterations nu=
	let (taillein,tailleout) = (vect_length in_data, vect_length out_data) in
	let erg=ref 1. and ergt = ref 1. and l=ref [] in
	if taillein<>tailleout then raise taille_incompatible;
	for iter = 0 to iterations-1 do
		ergt:=0.;
		for train = 0 to taillein-1 do
			entree N in_data.(train);
			propagation (N,P);
			erg:=erreurglobale N out_data.(train);
			modifpoids (N,P) (creergradient (N,P) out_data.(train)) (!erg) nu;
			ergt:=!ergt+. !erg;
		done;
		if (iter mod resolution) = 0 then l:=!ergt::(!l);
	done;
	rev !l;
;;

let apprentissage2 (N,P) (in_data,out_data) marge =
	let (taillein,tailleout) = (vect_length in_data, vect_length out_data) in
	if taillein<>tailleout then raise taille_incompatible;
	let erg = ref 1. and ergt = ref 1. and l = ref [] in
	while !ergt > marge do
		ergt:=0.;
		for train = 0 to taillein-1 do
			entree N in_data.(train);
			propagation (N,P);
			erg := erreurglobale N out_data.(train);
			modifpoids (N,P) (creergradient (N,P) out_data.(train)) !erg nu;
			ergt:=!ergt+. !erg;
		done;
	done;
;;

let showus (N,P) (in_data,out_data) =
	let (taillein,tailleout) = (vect_length in_data, vect_length out_data) in
	if taillein<>tailleout then raise taille_incompatible;
	for train = 0 to taillein-1 do
		print_floatvect in_data.(train);
		print_string " renvoie ";
		entree N in_data.(train);
		propagation (N,P);
		print_floatvect N.(tailleH-1);
		print_string " pour un résultat attendu de ";
		print_floatvect out_data.(train);
		print_string " soit une erreur globale de ";
		print_float (erreurglobale N out_data.(train));
		print_newline ();
	done;
;;

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
    trace_dot ((x2,y2)::q) mult;
;;

let set_or = ([|[|-1.;1.;1.|];[|-1.;0.;0.|];[|-1.;0.;1.|];[|-1.;1.;0.|]|],[|[|-1.;1.|];[|-1.;0.|];[|-1.;1.|];[|-1.;1.|]|]);;

let demontre sett tabnu =
	clear_graph();
	trace_cadre 1000 600 2 30 false;
	let n = vect_length tabnu in
	for i = 0 to n-1 do
		set_color (0xff5500/(i+1) + i*0x000011);
		let (N,P) = initpmc () in
      let bs = apprentissage1 (N,P) sett (1000*resolution) tabnu.(i) in
      let dots = conversion bs 1 in
      trace_dot dots 600.;
    done;
;;

demontre set_or [|0.1;0.15;0.20;0.25;0.3|];; (*EXCELLENT !!!!!!*)
