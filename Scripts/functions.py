MA	       = 'v($) = VM($) * MAR_RATE(S,NS,P,NP,Keq($))'
MA_deriv   = 'dvds = MAR_RATE_deriv(VM($),S,NS,P,NP,Keq($),&)'
MAIR       = 'v($) = VM($) * MAIR_RATE(S,NS)'
MAIR_deriv = 'dvds = MAIR_RATE_deriv(VM($),S,NS,&)'
MM	       = 'v($) = VM($) * MMR_RATE(S,NS,KS,P,NP, KP, Keq($))'
MM_deriv   = 'dvds = MMR_RATE_deriv(VM($),S,NS,KS,P,NP,KP,Keq($),&)'
MMIR       = 'v($) = VM($) * MMIR_RATE(S,NS,KS)'
MMIR_deriv = 'dvds = MMIR_RATE_deriv(VM($),S,NS,KS,&)'

Inhibit    = 'dfdi = RATE_INHIBIT_deriv(I,KIX,nh,&)'
Inhibit2   = 'dvds + v($)*dfdi;\r\n'
Activate   = 'dfda = RATE_ACTIVATE_deriv(A,KAX,nh,&)'
Activate2  = 'dvds + v($)*dfdi;\r\n'

Rate_Activate       = 'REGA = RATE_ACTIVATE(A,KAX,nh);\r\n'
Rate_Activate_deriv = 'RATE_ACTIVATE_deriv(A,KAX,nh,tow);\r\n'
Rate_Inhibit        = 'REGI = RATE_INHIBIT(I,KIX,nh);\r\n'
Rate_Inhibit_deriv  = 'RATE_INHIBIT_deriv(I,KIX,nh,tow);\r\n'
Rate_Act_Inh        = 'REGI = RATE_INHIBIT(I,KIX,nh);\r\n'+\
				      'REGA = RATE_ACTIVATE(A,KAX,nh);\r\n'
Rate_Act_Inh_deriv  = 'REGI = RATE_INHIBIT_deriv(I,KIX,nh);\r\n'+\
				      'REGA = RATE_ACTIVATE_deriv(A,KAX,nh);\r\n'
