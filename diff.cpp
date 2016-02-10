37d36
< 
39,40c38
<  // Nnumerov = 200;
<   Nnumerov = 400;
---
>   Nnumerov = 200;
44d41
< 
107,108c104,105
< //  Kwave2 = kconstant*pow(A/(1.+ A + Ecm/m0),2)*
< //		       Elab*(Elab/2./m0 + 1.);  
---
>   Kwave2 = kconstant*pow(A/(1.+ A + Ecm/m0),2)*
> 		       Elab*(Elab/2./m0 + 1.);  
110c107
<   Kwave2 = kconstant*mu*Ecm;
---
>   //Kwave2 = kconstant*mu*Ecm;
113,117d109
<   
< //cout<<"E =" << Ecm<<"\t Kwave in scatterRspace is" <<Kwave<<endl;
< //cout<<"kconstant =" << kconstant << "mu is " << mu<<endl;
< //cout<< "Mass"<<kconstant*mu<<endl;
<   gammaRel = 1;//2.*(Ecm/m0 +1.)/(Ecm/m0 + 2.);
118a111
>   gammaRel = 2.*(Ecm/m0 +1.)/(Ecm/m0 + 2.);
120d112
< //  cout<<"muhbar "<<muhbar<<endl;
137d128
<   
138a130
> 
139a132
>   
161c154
<            
---
> 
231c224,228
<       // 	cout<<Ecm<<" "<< "Elab "<<Elab << "  "<<l <<" "<<j<<" "<<SReal<<"  "<<SImag <<endl;
---
> 	  /*
> 	  if (Ecm > 190.) cout << l << " " << j << " " << S[l][up] << " " 
> 			       << 1. - S2[l][up] << " " << gg
>                                << " " << ff << endl;
> 	  */
233,243d229
<    	//  if (Ecm >50. )
<    	  //   {
<           //  	cout<<Ecm<<" "<<  l <<" "<<j<<" "<<SReal<<"  "<<SImag <<endl;
<           //  } 
<    	//  if (l==0 & j==.5 )
<    	    // {
<    	   if (Ecm>19 && Ecm <21 ) 
<                cout<<"Elab = " << Elab<< " Energ = " <<Ecm<<" L = " << l 
<                    <<" J =  "<<j<<" S real =  "<<SReal<<" S imag =  "<<SImag
<                   <<" S^2  "<<S2[l][up]<<" "<<" phase shift up = "<< phaseShift[l][up]<<endl;
<             //} 
291,292c277,278
<   //  for (int itry = 1;itry<20;itry++)
<     for (int itry = 1;itry<15;itry++)
---
> 
>     for (int itry = 1;itry<4;itry++)
295a282,285
> 	
> 	
> 
> 	
315,316d304
<  // complex<double> U (-Pot->HartreeFock.Vvol-Pot->VolumeAbove.Vvol,
<    //                                        -Pot->VolumeAbove.Wvol);
321,322c309
< //  for (int i=0;i<1000;i++)
<   for (int i=0;i<1000;i++)
---
>   for (int i=0;i<10;i++)
506,507c493,494
<    
<    
---
> 
> 
522d508
< 
525d510
< 
527,528d511
< 
< 
540d522
< 
