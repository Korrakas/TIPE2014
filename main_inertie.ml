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
	let neur, poids, souv = make_vect tailleH [||], make_vect tailleH [||], make_vect tailleH [||]
	in
	let fillp pij t =
		for k = 1 to t-1 do (*pas de lien vers le biais*)
			pij.(k) <- random__float 0.1 -. 0.05; (* on remplit le tableau des poids synaptiques aléatoirement entre 0.5 et 1.*)
		done;
	in
	for i = 0 to tailleH - 1 do
		neur.(i) <- make_vect config.(i) 0.;
		poids.(i) <- make_vect config.(i) [||];
		souv.(i) <- make_vect config.(i) [||];
		neur.(i).(0) <- -1.;
		if i <> tailleH - 1 then (*pas de liens externes pour la couche de sortie*) 
		for j = 0 to config.(i) - 1 do
			poids.(i).(j) <- make_vect config.(i+1) 0.;
			souv.(i).(j) <- make_vect config.(i+1) 0.;
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

let nu = 0.1;;

let alpha = 0.;;

let modifpoids_inertiel (neur,poids,souvenir) grad erg =
	for i = 0 to tailleH - 2 do
		for j = 0 to config.(i) - 1 do
			for h = 1 to config.(i+1) -1 do (*pas besoin de toucher le "poids vers le biais" qui ne doit de toute façon pas exister (0 ?) *)
				souvenir.(i).(j).(h) <- nu*.grad.(i+1).(h)*.(actprime erg)*.neur.(i).(j) + alpha*.souvenir.(i).(j).(h);
				poids.(i).(j).(h) <- poids.(i).(j).(h) +. souvenir.(i).(j).(h); 
			done;
		done;
	done;
	memoire;
;;

let apprentissage1_inertiel (N,P,S) (in_data,out_data) iterations =
	let (taillein,tailleout) = (vect_length in_data, vect_length out_data) in
	let erg=ref 1. and ergt = ref 1. and l=ref [] in
	let souv = ref (map_vect copy_vect P) in
	if taillein<>tailleout then raise taille_incompatible;
	for iter = 0 to iterations-1 do
		ergt:=0.;
		for train = 0 to taillein-1 do
			entree N in_data.(train);
			propagation (N,P);
			erg:=erreurglobale N out_data.(train);
			modifpoids_inertiel (N,P,S) (creergradient (N,P) out_data.(train)) (!erg);
			ergt:=!ergt+. !erg;
		done;
		l:=!ergt::(!l);
	done;
	rev !l;
;;

let apprentissage2_inertiel (N,P,S) (in_data,out_data) marge =
	let (taillein,tailleout) = (vect_length in_data, vect_length out_data) in
	if taillein<>tailleout then raise taille_incompatible;
	let erg = ref 1. and ergt = ref 1. and l = ref [] in
	while !ergt > marge do
		ergt:=0.;
		for train = 0 to taillein-1 do
			entree N in_data.(train);
			propagation (N,P);
			erg := erreurglobale N out_data.(train);
			modifpoids (N,P) (creergradient (N,P) out_data.(train)) !erg;
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


(* Routine test OR*)
let (N,P,S) = initpmc ();;
showus (N,P) ([|[|-1.;1.;1.|];[|-1.;0.;0.|];[|-1.;0.;1.|];[|-1.;1.;0.|]|],[|[|-1.;1.|];[|-1.;0.|];[|-1.;1.|];[|-1.;1.|]|]);;
apprentissage2 (N,P,S) ([|[|-1.;1.;1.|];[|-1.;0.;0.|];[|-1.;0.;1.|];[|-1.;1.;0.|]|],[|[|-1.;1.|];[|-1.;0.|];[|-1.;1.|];[|-1.;1.|]|]) 0.001;;
showus (N,P) ([|[|-1.;1.;1.|];[|-1.;0.;0.|];[|-1.;0.;1.|];[|-1.;1.;0.|]|],[|[|-1.;1.|];[|-1.;0.|];[|-1.;1.|];[|-1.;1.|]|]);;

(* Routine test NAND *)
let (N,P,S) = initpmc ();;
showus (N,P) ([|[|-1.;1.;1.|];[|-1.;0.;0.|];[|-1.;0.;1.|];[|-1.;1.;0.|]|],[|[|-1.;0.|];[|-1.;1.|];[|-1.;1.|];[|-1.;1.|]|]);;
apprentissage2 (N,P,S) ([|[|-1.;1.;1.|];[|-1.;0.;0.|];[|-1.;0.;1.|];[|-1.;1.;0.|]|],[|[|-1.;0.|];[|-1.;1.|];[|-1.;1.|];[|-1.;1.|]|]) 0.001;;
showus (N,P) ([|[|-1.;1.;1.|];[|-1.;0.;0.|];[|-1.;0.;1.|];[|-1.;1.;0.|]|],[|[|-1.;0.|];[|-1.;1.|];[|-1.;1.|];[|-1.;1.|]|]);;
