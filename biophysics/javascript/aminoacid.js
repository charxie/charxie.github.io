//JavaScript Slideshow for 20 Animo Acids
//(C) Qian Xie, 1999, University of Cyprus

 var index = 1;
 var selected = 1;
 var total = 20;
 
 caption=new Array(total);
 caption[ 1] = "Alanine" 
 caption[ 2] = "Arginine" 
 caption[ 3] = "Asparagine" 
 caption[ 4] = "Aspartate" 
 caption[ 5] = "Cysteine" 
 caption[ 6] = "Glutamine" 
 caption[ 7] = "Glutamate" 
 caption[ 8] = "Glycine" 
 caption[ 9] = "Histidine" 
 caption[10] = "Isoleucine" 
 caption[11] = "Leucine" 
 caption[12] = "Lysine" 
 caption[13] = "Methionine" 
 caption[14] = "Phenylalanine" 
 caption[15] = "Proline" 
 caption[16] = "Serine" 
 caption[17] = "Threonine" 
 caption[18] = "Tryptophan" 
 caption[19] = "Tyrosine" 
 caption[20] = "Valine" 

 oneLetter=new Array(total);
 oneLetter[ 1]="A";
 oneLetter[ 2]="R";
 oneLetter[ 3]="N";
 oneLetter[ 4]="D";
 oneLetter[ 5]="C";
 oneLetter[ 6]="Q";
 oneLetter[ 7]="E";
 oneLetter[ 8]="G";
 oneLetter[ 9]="H";
 oneLetter[10]="I";
 oneLetter[11]="L";
 oneLetter[12]="K";
 oneLetter[13]="M";
 oneLetter[14]="F";
 oneLetter[15]="P";
 oneLetter[16]="S";
 oneLetter[17]="T";
 oneLetter[18]="W";
 oneLetter[19]="Y";
 oneLetter[20]="V";

 charge=new Array(total);
 charge[ 1]="0";
 charge[ 2]="+";
 charge[ 3]="0";
 charge[ 4]="-";
 charge[ 5]="0";
 charge[ 6]="0";
 charge[ 7]="-";
 charge[ 8]="0";
 charge[ 9]="+";
 charge[10]="0";
 charge[11]="0";
 charge[12]="+";
 charge[13]="0";
 charge[14]="0";
 charge[15]="0";
 charge[16]="0";
 charge[17]="0";
 charge[18]="0";
 charge[19]="0";
 charge[20]="0";

 sideChain=new Array(total);
 sideChain[ 1]="aliphatic";
 sideChain[ 2]="basic";
 sideChain[ 3]="amide";
 sideChain[ 4]="acidic";
 sideChain[ 5]="surfur-containing";
 sideChain[ 6]="amide";
 sideChain[ 7]="acidic";
 sideChain[ 8]="aliphatic";
 sideChain[ 9]="basic";
 sideChain[10]="aliphatic";
 sideChain[11]="aliphatic";
 sideChain[12]="basic";
 sideChain[13]="surfur-containing";
 sideChain[14]="aromatic";
 sideChain[15]="secondary amino group";
 sideChain[16]="aliphatic hydroxyl";
 sideChain[17]="aliphatic hydroxyl";
 sideChain[18]="aromatic";
 sideChain[19]="aromatic";
 sideChain[20]="aliphatic";

 function getPicture(selected){
   document.aminoacid.src="aminoacid/"+caption[selected]+".gif";
   document.acids.caption.value = caption[selected];
   document.acids.oneletter.value = oneLetter[selected];
   document.acids.charge.value = charge[selected];
   document.acids.sidechain.value = sideChain[selected];
   index=selected;
 }

 function decrement(form){ 
   index --;
   if (index == 0){ 
   index = total;
   } 
   document.aminoacid.src="aminoacid/"+caption[index]+".gif";
   form.caption.value = caption[index];
   document.acids.oneletter.value = oneLetter[index];
   document.acids.charge.value = charge[index];
   document.acids.sidechain.value = sideChain[index];
 
 } 
 function increment(form){ 
   index ++;
   if (index == total+1 ){ 
   index = 1; 
   }
   document.aminoacid.src="aminoacid/"+caption[index]+".gif";
   form.caption.value = caption[index];
   document.acids.oneletter.value = oneLetter[index];
   document.acids.charge.value = charge[index];
   document.acids.sidechain.value = sideChain[index];
 }
