function v = RATE_INHIBIT(I,KIX,nh); % a factor for inhibtion 
    II = 1 + (I./KIX).^nh; 
    v = prod(1./II); 
end;

function dfdi = RATE_INHIBIT_deriv(I,KIX,nh,tow); % a factor for inhibtion 
% the function takes a vector of inhibitors and calculates the product of hill functions
    II = 1 + (I./KIX).^nh; 
    FA = prod(1./II); 
	ia = (KIX(tow)/I(tow))^nh(tow);
	fact = 1./(1+ia);
	dfdi = -1*nh(tow)*FA*fact/I(tow); 
end;

function MA = MAIR_RATE(S,NS); % a reversible rate equation
    SN = S.^NS;
    MA = prod(SN);
end;

% the function provides derivative for a irreversible
% mass-action equation
% Input are substrate concentrations, Keq
% tow: index of derivative 
% WARNING: THE FUNCTION MAY HAVE PROBLEMS WHEN RATE OR
% CONCENTRATION IS ZERO .... CHECK!

function dvds = MAIR_RATE_deriv(VM,S,NS,tow); % a reversible rate equation
	SN = S.^NS;
	v = VM*prod(SN);

	dvds = NS(tow)*v/S(tow);
end;

function v = MAR_RATE(S,NS,P,NP,Keq) % a reversible rate equation;
    SN = S.^NS; PN = P.^NP;
    v = (prod(SN) - prod(PN)/Keq);
end;

% the function provides derivative for a reversible
% mass-action equation
% Input is substrate/product concentrations, Keq
% tow: index of derivative (negative for products)
% WARNING: THE FUNCTION MAY HAVE PROBLEMS WHEN RATE OR
% CONCENTRATION IS ZERO .... CHECK!

function dvds = MAR_RATE_deriv(VM,S,NS,P,NP,Keq,tow); % a reversible rate equation
	SN = S.^NS; PN = P.^NP;
	vforward = VM*prod(SN);
	vback = VM*prod(PN)/Keq;

	if (tow>0); 
		dvds = NS(tow)*vforward/S(tow);
	end;

	if (tow<0);
		tow=abs(tow);
		dvds = -1*NP(tow)*vback/P(tow);
	end;

	if (tow==0), error("derivative incorrect"), end;
end;

function v = MMIR_RATE(S,NS,KS); % a irreversible rate equation
    MA = prod((S./KS).^NS);
    SX = prod((1+S./KS).^NS);
    v  = MA/SX;
end;

% the function provides derivative for a generic irreversible
% Michaelis-Menten equation
% tow: index of derivative 
% WARNING: THE FUNCTION MAY HAVE PROBLEMS WHEN RATE OR
% CONCENTRATION IS ZERO .... CHECK!

function dvds = MMIR_RATE_deriv(VM,S,NS,KS,tow); % a irreversible rate equation
	MA = prod((S./KS).^NS);
	SX = prod((1+S./KS).^NS);
	v0  = VM*MA/SX;

	deri = 1/(1+S(tow)/KS(tow));
	dvds = NS(tow)*v0*deri/S(tow);
end;

function v = MMR_RATE(S,NS,KS,P,NP,KP,Keq); % a reversible rate equation
    SN = S.^NS; PN = P.^NP;
    MA = (prod(SN) - prod(PN)/Keq);

    SX = (1+S./KS).^NS; PX = (1+ P./KP).^NP;
    F = prod(SX) + prod(PX) - 1;
    KX = 1/prod(KS.^NS);

    v = KX*MA/F;
end;

% the function provides the rate equation for a generic reversible Michaelis-Menten
% equation
% Input is substrate/product concentrations, Keq, and kinetic coefficients

function dvds = MMR_RATE_deriv(VM,S,NS,KS,P,NP,KP,Keq,tow); % a reversible rate equation

    SN = S.^NS; PN = P.^NP;
    MA = (prod(SN) - prod(PN)/Keq);

    SX = (1+S./KS).^NS;
    PX = (1+ P./KP).^NP;
    F = prod(SX) + prod(PX) - 1;
    KX = 1/prod(KS.^NS);

    v0 = VM*KX*MA/F;
	
	v_plus  = VM*KX*prod(SN)/F;
	v_minus = VM*KX*prod(PN)/(F*Keq);
	ratio = (1 - prod(PN)/(prod(SN)*Keq));

	% the derivative
	if (tow>0);
		deri = (S(tow)/KS(tow))/(1 + S(tow)/KS(tow));
		dvds = (NS(tow)*v_plus/S(tow))*(1 - prod(SX) * deri * ratio / F);
	end;
	
	if (tow<0);
		tow = abs(tow);
		ratio = (v_plus/v_minus - 1);
		ppx = P(tow)/KP(tow);
		deri = ppx/(1+ppx);
		dvds = -1*(NP(tow)*v_minus/P(tow))*(1 + prod(PX) * deri * ratio / F);
	end;
end;

function v = RATE_ACTIVATE(A,KAX,nh); 
    AA = (A./KAX).^nh; 
    v = prod(AA./(1+AA)); 
end;

function dfda = RATE_ACTIVATE_deriv(A,KAX,nh,tow); % a pre-factor for activation
% the function works on vectors 
    AA = (A./KAX).^nh; 
    FA = prod(AA./(1+AA)); 

	fact = 1 / ((A(tow)/KAX(tow))^nh(tow));
	dfda =  nh(tow)*FA*fact/A(tow);
end;
