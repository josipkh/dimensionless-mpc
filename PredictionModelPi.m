function out1 = PredictionModelPi(Ts,in2)
%PREDICTIONMODELPI
%    OUT1 = PREDICTIONMODELPI(TS,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    28-Jul-2021 15:00:34

Tflp = in2(8,:);
Tfrp = in2(9,:);
Trlp = in2(10,:);
Trrp = in2(11,:);
deltaf = in2(12,:);
thetadp = in2(3,:);
vxp = in2(1,:);
vyp = in2(2,:);
wflp = in2(4,:);
wfrp = in2(5,:);
wrlp = in2(6,:);
wrrp = in2(7,:);
t2 = deltaf./1.3e+1;
t5 = sqrt(8.74e+2);
t6 = 5.860887304837041e+1;
t20 = vxp.*2.079219468935398e+5;
t22 = vxp.*1.299512168084624e+1;
t23 = vyp.*1.299512168084624e+1;
t24 = vxp.*1.299512168084624e+4;
t26 = vyp.*1.299512168084624e+1i;
t27 = vyp.*1.299512168084624e+4i;
t3 = cos(t2);
t4 = sin(t2);
t10 = t5.*t6.*thetadp.*6.0e+1;
t11 = t5.*t6.*thetadp.*(3.0./8.0e+2);
t13 = t5.*t6.*thetadp.*(3.75+3.75i);
t14 = t5.*t6.*thetadp.*(-3.75+3.75i);
t15 = t5.*t6.*thetadp.*(3.75-3.75i);
t19 = t5.*t6.*thetadp.*(3.75e-3+3.75e-3i);
t21 = t5.*t6.*thetadp.*(3.75e-3-3.75e-3i);
t7 = t4.*1i;
t12 = -t11;
t16 = t3.*t5.*t6.*thetadp.*(1.5e+1./4.0);
t17 = t4.*t5.*t6.*thetadp.*(1.5e+1./4.0);
t25 = -t19;
t28 = t3.*t24;
t29 = t4.*vyp.*1.299512168084624e+4;
t30 = t11+t22;
t31 = t11+t23;
t34 = -1.0./(t11-t22);
t37 = -t3.*(t11-t22);
t42 = t4.*(t11-t22).*1i;
t43 = t13+t24+t27;
t44 = t14+t24+t27;
t45 = t21+t22+t26;
t8 = -t7;
t18 = -t16;
t32 = t12+t22;
t33 = 1.0./t30;
t35 = t3.*t30;
t36 = t4.*t31;
t38 = t3.*t31.*1i;
t39 = t7.*t30;
t40 = -t7.*(t11-t22);
t41 = t4.*t30.*-1i;
t46 = angle(t45);
t47 = t22+t25+t26;
t57 = t16+t17+t28+t29;
t9 = t3+t8;
t48 = angle(t47);
t53 = t35+t36;
t54 = t36+t37;
t58 = t17+t18+t28+t29;
t59 = 1.0./t57;
t49 = t9.*t43;
t51 = t9.*t44;
t55 = 1.0./t53;
t56 = 1.0./t54;
t60 = 1.0./t58;
t65 = t38+t41+t53;
t67 = t38+t42+t54;
t50 = angle(t49);
t52 = angle(t51);
t61 = t5.*t6.*t55.*wfrp.*9.37354118993135e+1;
t62 = t5.*t6.*t56.*wflp.*9.37354118993135e+1;
t66 = angle(t65);
t68 = angle(t67);
t63 = t61-1.0305e+5;
t64 = t62-1.0305e+5;
out1 = [t5.*t6.*(t3.*(-1.288125e+2)+t4.*t50.*4.1901875e+1+t4.*t52.*4.1901875e+1+t5.*t6.*thetadp.*vyp.*3.717140068891945e-2+(t5.*t6.*wrrp.*9.37354118993135e+2)./(t10+t20)-(t5.*t6.*wrlp.*9.37354118993135e+2)./(t10-t20)+t3.*t5.*t6.*t60.*wflp.*5.858463243707094e+1+t3.*t5.*t6.*t59.*wfrp.*5.858463243707094e+1-1.288125e+2).*8.960926177932477e-6;t5.*t6.*(t4.*1.288125e+2+t46.*4.1901875e+1+t48.*4.1901875e+1+t3.*t50.*4.1901875e+1+t3.*t52.*4.1901875e+1+t5.*t6.*thetadp.*vxp.*3.717140068891945e-2-t4.*t5.*t6.*t60.*wflp.*5.858463243707094e+1-t4.*t5.*t6.*t59.*wfrp.*5.858463243707094e+1).*(-8.960926177932477e-6);t46.*1.494849379741506+t48.*1.494849379741506+t3.*t63.*4.459374967532795e-5-t3.*t64.*4.459374967532795e-5+t4.*t63.*2.229687483766398e-5+t4.*t64.*2.229687483766398e-5-t3.*t66.*1.494849379741506+t4.*t66.*2.989698759483012-t3.*t68.*1.494849379741506-t4.*t68.*2.989698759483012+t5.*t6.*t33.*wrrp.*4.180013493951743e-3+(t5.*t6.*wrlp.*4.180013493951743e-3)./(t11-t22);Tflp.*9.016241311475409e+3-t5.*t6.*t56.*wflp.*9.946622950819672e-1+1.093502950819672e+3;Tfrp.*9.016241311475409e+3-t5.*t6.*t55.*wfrp.*9.946622950819672e-1+1.093502950819672e+3;Trlp.*9.016241311475409e+3+(t5.*t6.*wrlp.*9.946622950819672e-1)./(t11-t22)+1.093502950819672e+3;Trrp.*9.016241311475409e+3-t5.*t6.*t33.*wrrp.*9.946622950819672e-1+1.093502950819672e+3;Tflp;Tfrp;Trlp;Trrp;deltaf];
