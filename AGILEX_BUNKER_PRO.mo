package AGILEX_BUNKER_PRO
  model bush_element
    Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a annotation (
      Placement(visible = true, transformation(origin = {-100, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-98, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b annotation (
      Placement(visible = true, transformation(origin = {98, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {100, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    parameter Real k = 100;
    parameter Modelica.Units.SI.ModulusOfElasticity Er = 1e7;
    parameter Modelica.Units.SI.PoissonNumber vs = 0.3;
    parameter Modelica.Units.SI.PoissonNumber vr = 0.5;
    parameter Modelica.Units.SI.Length h = 0.012;
    parameter Modelica.Units.SI.Length L = 0.0586*38;
    parameter Modelica.Units.SI.Length b = 0.15;
    final parameter Modelica.Units.SI.Length l = L/80;
    final parameter Modelica.Units.SI.Area A = h*b;
    final parameter Modelica.Units.SI.ShearModulus Gr = Er/(2*(1+vr));
    final parameter Modelica.Units.SI.TranslationalSpringConstant Kdx = ((Er*A)/l)*k;//1e7;
    final parameter Modelica.Units.SI.TranslationalDampingConstant rdx = Kdx*0.1;//1e6;
    final parameter Modelica.Units.SI.TranslationalSpringConstant Kdy = ((Gr*A)/l)*k;//1e7;
    final parameter Modelica.Units.SI.TranslationalDampingConstant rdy = Kdy*0.1;//1e6;
    final parameter Modelica.Units.SI.TranslationalSpringConstant Kdz = ((Gr*A)/l)*k;//1e7;
    final parameter Modelica.Units.SI.TranslationalDampingConstant rdz = Kdz*0.1;//1e6;
    //10 ^ 5;
    final parameter Modelica.Units.SI.RotationalSpringConstant Kx = ((Gr*h*b^3)/(12*l))*k;//1e4;
    //573;
    final parameter Modelica.Units.SI.RotationalSpringConstant Ky = ((Er*h*b^3)/(12*l))*k;//1e4;
    final parameter Modelica.Units.SI.RotationalSpringConstant Kz = ((Er*b*h^3)/(12*l));// + ((Es*b*h^3)/(12*l))*0.1^3;//10;
    //11.46;
    final parameter Modelica.Units.SI.RotationalDampingConstant rx = Kx*0.1;//1e3;
    //5.73;
    final parameter Modelica.Units.SI.RotationalDampingConstant ry = Ky*0.1;//1e3;
    final parameter Modelica.Units.SI.RotationalDampingConstant rz = Kz*0.01;//0.1;
    //0.1146;
    Modelica.Units.SI.Angle rth[3], rth1[3];
    Modelica.Mechanics.MultiBody.Frames.Orientation Rrel;
    Real E[3,3];
  equation
    //Forza Elastica
    frame_a.f = -{Kdx, Kdy, Kdz} .* Modelica.Mechanics.MultiBody.Frames.resolve2(frame_a.R, (frame_b.r_0 - frame_a.r_0)) - {rdx, rdy, rdz} .* Modelica.Mechanics.MultiBody.Frames.resolve2(frame_a.R, (der(frame_b.r_0) - der(frame_a.r_0)));
    //Equilibrio delle forze
    frame_b.f + Modelica.Mechanics.MultiBody.Frames.resolve2(Rrel, frame_a.f) = {0, 0, 0};
    //Rotazione relativa
    Rrel = Modelica.Mechanics.MultiBody.Frames.relativeRotation(frame_a.R, frame_b.R);
    E = [1, (Modelica.Math.sin(rth[1])*Modelica.Math.sin(rth[2]))/Modelica.Math.cos(rth[2]), -(Modelica.Math.cos(rth[1])*Modelica.Math.sin(rth[2]))/Modelica.Math.cos(rth[2]);
       0, Modelica.Math.cos(rth[1]), Modelica.Math.sin(rth[1]);
       0, -Modelica.Math.sin(rth[1])/Modelica.Math.cos(rth[2]), Modelica.Math.cos(rth[1])/Modelica.Math.cos(rth[2])];

    //Calcolo dell'angolo relativo rispetto al frame_a
    der(rth) = E*Modelica.Mechanics.MultiBody.Frames.resolve1(Rrel, Rrel.w);

    rth1 = Modelica.Mechanics.MultiBody.Frames.axesRotationsAngles(Rrel, {1,2,3});

    when Modelica.Math.Vectors.length(rth - rth1) > 1e-5 then
    reinit(rth, rth1);
    end when;

    //Torsione elestica
    frame_a.t = - transpose(E)*{Kx * rth[1] + rx * der(rth[1]), Ky * rth[2] + ry * der(rth[2]), Kz * rth[3] + rz * der(rth[3])};
    //Equilibrio torsionale
    frame_b.t + Modelica.Mechanics.MultiBody.Frames.resolve2(Rrel, frame_a.t) = {0, 0, 0};
    // + cross(Modelica.Mechanics.MultiBody.Frames.resolve2(frame_b.R, (frame_a.r_0 - frame_b.r_0)), Modelica.Mechanics.MultiBody.Frames.resolve2(Rrel, frame_a.f))
    //Assegnazione angolo iniziale
  initial equation
    rth =  Modelica.Mechanics.MultiBody.Frames.axesRotationsAngles(Rrel, {1,2,3});
    //rth = {Modelica.Math.atan2(-Rrel.T[3, 2], Rrel.T[3, 3]), Modelica.Math.asin(Rrel.T[3, 1]), Modelica.Math.atan2(-Rrel.T[2, 1], Rrel.T[1, 1])}; //usarlo quando si avvia track_bunker_pro e track_bunker_pro_new
    annotation (
      Icon(graphics={  Ellipse(fillColor = {255, 0, 0}, fillPattern = FillPattern.Solid, extent = {{-100, 100}, {100, -100}})}));
  end bush_element;

  model element
    parameter Modelica.Units.SI.Length startX = 0;
    parameter Modelica.Units.SI.Length startY = 0;
    parameter Modelica.Units.SI.Length startZ = 0;
    parameter Modelica.Units.SI.Angle thX = 0;
    parameter Modelica.Units.SI.Angle thY = 0;
    parameter Modelica.Units.SI.Angle thZ = 0;
    parameter Real l = 2.28;
    Modelica.Mechanics.MultiBody.Parts.BodyShape bodyShape(angles_fixed = true, angles_start(each displayUnit = "rad") = {thX, thY, thZ}, animateSphere = false, height = 0.15,
      length=l/80,
      m=2.5/80,
      r={l/80,0,0},                                                                                                                                                                                                        r_0(each fixed = true, start = {startX, startY, startZ}),
      r_CM={(l/80)*0.5,0,0},                                                                                                                                                                                                        shapeType = "box", useQuaternions = false,
      width=0.0142)                                                                                                                                                                                                         annotation (
      Placement(visible = true, transformation(origin = {-12, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a annotation (
      Placement(visible = true, transformation(origin = {-92, 0}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-92, 0}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a1 annotation (
      Placement(visible = true, transformation(origin = {-8, -74}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-4, -98}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation(animation = false, r={0,-0.0071,
          -0.15*0.5})                                                                                                      annotation (
      Placement(visible = true, transformation(origin = {-10, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation1(animation = false, r={(l/80)
          *0.5,0.0071,0})                                                                                                   annotation (
      Placement(visible = true, transformation(origin = {-26, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    AGILEX_BUNKER_PRO.bush_element bush_element annotation (Placement(visible=
            true, transformation(
          origin={38,4},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b1 annotation (
      Placement(visible = true, transformation(origin = {-12, 72}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(extent = {{-20, 76}, {12, 108}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b annotation (
      Placement(visible = true, transformation(origin = {90, 6}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Units.SI.Velocity v0[3];
    Modelica.Units.SI.Position r0[3];
    Modelica.Units.SI.Angle th0;
  equation
    v0 = der(bodyShape.frame_a.r_0);
    r0 = bodyShape.frame_a.r_0;
    th0 = Modelica.Math.atan2(bodyShape.frame_b.r_0[2] - bodyShape.frame_a.r_0[2],bodyShape.frame_b.r_0[1] - bodyShape.frame_a.r_0[1]);
    connect(frame_a, bodyShape.frame_a) annotation (
      Line(points = {{-92, 0}, {-22, 0}, {-22, 4}}));
    connect(bodyShape.frame_a, fixedTranslation.frame_a) annotation (
      Line(points = {{-22, 4}, {-36, 4}, {-36, -32}, {-20, -32}}, color = {95, 95, 95}));
    connect(fixedTranslation1.frame_a, bodyShape.frame_a) annotation (
      Line(points = {{-36, 40}, {-56, 40}, {-56, 4}, {-22, 4}}, color = {95, 95, 95}));
    connect(bodyShape.frame_b, bush_element.frame_a) annotation (
      Line(points={{-2,4},{14,4},{14,3.8},{28.2,3.8}}));
    connect(fixedTranslation1.frame_b, frame_b1) annotation (
      Line(points = {{-16, 40}, {-12, 40}, {-12, 72}}));
    connect(bush_element.frame_b, frame_b) annotation (
      Line(points={{48,3.8},{90,3.8},{90,6}},    color = {95, 95, 95}));
    connect(fixedTranslation.frame_b, frame_a1) annotation (
      Line(points = {{0, -32}, {76, -32}, {76, -74}, {-8, -74}}));
    annotation (
      Icon(graphics={  Rectangle(origin = {-6, -1}, fillColor = {113, 113, 113}, fillPattern = FillPattern.Solid, extent = {{-34, 99}, {34, -99}})}),
      Diagram);
  end element;

  model track
    extends Modelica.Icons.Example;
    inner Modelica.Mechanics.MultiBody.World world annotation (
      Placement(visible = true, transformation(origin = {-368, 186}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder(diameter = 0.16, length = 0.1, lengthDirection = {0, 0, 1}, r = {0, 0, 0.05}) annotation (
      Placement(visible = true, transformation(origin = {-120, 276}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute(animation = false, useAxisFlange = true) annotation (
      Placement(visible = true, transformation(origin = {-186, 270}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder4(diameter = 0.088, length = 0.075, r = {0, 0, 0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-124, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute1 annotation (
      Placement(visible = true, transformation(origin = {-167, 55}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder1(diameter = 0.088, length = 0.075, r = {0, 0, -0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-124, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder2(diameter = 0.088, length = 0.075, r = {0, 0, 0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-126, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute2 annotation (
      Placement(visible = true, transformation(origin = {-157, 11}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder5(diameter = 0.088, length = 0.075, r = {0, 0, -0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-129, -35}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder6(diameter = 0.088, length = 0.075, r = {0, 0, 0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-132, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute3 annotation (
      Placement(visible = true, transformation(origin = {-173, -67}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder7(diameter = 0.130, length = 0.075, r = {0, 0, 0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-124, 114}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder8(diameter = 0.130, length = 0.075, r = {0, 0, -0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-126, 158}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute4 annotation (
      Placement(visible = true, transformation(origin = {-169, 127}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.Fixed fixed(animation = false, r = {0.858, 0.065, 0.05}) annotation (
      Placement(visible = true, transformation(origin = {-234, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation(animation = true, r = {-0.153, 0.040, 0}) annotation (
      Placement(visible = true, transformation(origin = {-198, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation1(r = {-0.350, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-188, 22}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation2(r = {-0.25, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-190, -6}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute5 annotation (
      Placement(visible = true, transformation(origin = {-261, -125}, extent = {{-13, 13}, {13, -13}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation3(animation = false, r = {0.250, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-300, -126}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation4(r = {-0.192, -0.192, 0}) annotation (
      Placement(visible = true, transformation(origin = {-218, -124}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder9(diameter = 0.105, length = 0.075, r = {0, 0, -0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-128, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute6 annotation (
      Placement(visible = true, transformation(origin = {-181, -121}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder10(diameter = 0.105, length = 0.075, r = {0, 0, 0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-128, -148}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder3(diameter = 0.088, length = 0.075, r = {0, 0, -0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-124, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder11(diameter = 0.105, length = 0.075, r = {0, 0, 0.0375}) annotation (
      Placement(visible = true, transformation(origin={-134,-226},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute7 annotation (
      Placement(visible = true, transformation(origin = {-185, -203}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation5(animation = false, r = {0.281, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-310, -222}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder12(diameter = 0.105, length = 0.075, r = {0, 0, -0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-134, -182}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute8 annotation (
      Placement(visible = true, transformation(origin = {-276, -224}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation7(r = {-0.087, -0.029, 0}) annotation (
      Placement(visible = true, transformation(origin = {-288, 212}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder13(diameter = 0.105, length = 0.075, r = {0, 0, 0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-122, 194}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute9 annotation (
      Placement(visible = true, transformation(origin = {-175, 217}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder14(diameter = 0.105, length = 0.075, r = {0, 0, -0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-122, 230}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation9(r = {-0.071, -0.120, 0}) annotation (
      Placement(visible = true, transformation(origin = {-208, 216}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute10 annotation (
      Placement(visible = true, transformation(origin = {-246, 218}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation6(r = {-0.192, -0.192, 0}) annotation (
      Placement(visible = true, transformation(origin = {-226, -202}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation8(animation = false, r = {-0.134, -0.134, 0}) annotation (
      Placement(visible = true, transformation(origin = {-288, -84}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation10(animation = false, r = {0.064, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-208, -48}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel springDamperParallel(
      c=100000,                                                                             d = 1000, fixedRotationAtFrame_a = false, fixedRotationAtFrame_b = false,
      s_unstretched=0.15)                                                                                                                                                                   annotation (
      Placement(visible = true, transformation(origin = {-264, -58}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation11(animation = false, r = {-0.134, -0.134, 0}) annotation (
      Placement(visible = true, transformation(origin = {-388, -182}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation12(animation = false, r = {0.098, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-374, -10}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel springDamperParallel1(
      c=100000,                                                                              d = 1000, fixedRotationAtFrame_a = false, fixedRotationAtFrame_b = false,
      s_unstretched=0.14)                                                                                                                                                                    annotation (
      Placement(visible = true, transformation(origin = {-401, -83}, extent = {{-11, -11}, {11, 11}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel springDamperParallel2(
      c=100000,                                                                              d = 1000, fixedRotationAtFrame_a = false, fixedRotationAtFrame_b = false,
      s_unstretched=0.17)                                                                                                                                                                    annotation (
      Placement(visible = true, transformation(origin = {-410, 56}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation13(animation = false, r = {-0.044, -0.065, 0}) annotation (
      Placement(visible = true, transformation(origin = {-250, 142}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation14(animation = false, r = {0.278, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-350, 30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    AGILEX_BUNKER_PRO.cinghia cinghia annotation (Placement(visible=true,
          transformation(
          origin={177,21},
          extent={{-41,41},{41,-41}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact_sprocket contact_sprocket annotation (Placement(
          visible=true, transformation(
          origin={28,240},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact6(R=0.0525) annotation (Placement(visible=
            true, transformation(
          origin={-142,-322},
          extent={{-32,-32},{32,32}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact5(R=0.0525) annotation (Placement(visible=
            true, transformation(
          origin={-51,-299},
          extent={{-33,-33},{33,33}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact4(R=0.044) annotation (Placement(visible=
            true, transformation(
          origin={-38,-210},
          extent={{-26,-26},{26,26}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact3(R=0.044) annotation (Placement(visible=
            true, transformation(
          origin={-19,-135},
          extent={{-27,-27},{27,27}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact2(R=0.044) annotation (Placement(visible=
            true, transformation(
          origin={-6,-44},
          extent={{-34,-34},{34,34}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact(R=0.0525) annotation (Placement(visible=
            true, transformation(
          origin={10,174},
          extent={{-30,-30},{30,30}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact1(R=0.065) annotation (Placement(visible=
            true, transformation(
          origin={12,68},
          extent={{-36,-36},{36,36}},
          rotation=0)));
    Modelica.Mechanics.Rotational.Sources.Position position annotation (
      Placement(visible = true, transformation(origin = {-412, 300}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.RealExpression realExpression(y = 0) annotation (
      Placement(visible = true, transformation(origin = {-528, 292}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Prismatic prismatic(useAxisFlange=true)
                       annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-210,94})));
    Modelica.Mechanics.Translational.Sources.Position position1(exact=true)
      annotation (Placement(transformation(extent={{-260,90},{-240,110}})));
    Modelica.Blocks.Sources.RealExpression realExpression1
      annotation (Placement(transformation(extent={{-314,58},{-294,78}})));
  equation
    connect(revolute.frame_a, world.frame_b) annotation (
      Line(points = {{-196, 270}, {-196, 264.25}, {-228, 264.25}, {-228, 236.5}, {-266, 236.5}, {-266, 241}, {-358, 241}, {-358, 186}}));
    connect(revolute.frame_b, bodyCylinder.frame_a) annotation (
      Line(points = {{-176, 270}, {-140, 270}, {-140, 276}, {-130, 276}}));
    connect(bodyCylinder4.frame_a, revolute1.frame_b) annotation (
      Line(points = {{-134, 60}, {-134, 52.75}, {-138, 52.75}, {-138, 55.5}, {-156, 55.5}, {-156, 55}}, color = {95, 95, 95}));
    connect(bodyCylinder1.frame_a, revolute2.frame_b) annotation (
      Line(points = {{-134, 28}, {-146, 28}, {-146, 11}}, color = {95, 95, 95}));
    connect(bodyCylinder2.frame_a, revolute2.frame_b) annotation (
      Line(points = {{-136, -4}, {-136, -7.5}, {-146, -7.5}, {-146, 11}}, color = {95, 95, 95}));
    connect(bodyCylinder6.frame_a, revolute3.frame_b) annotation (
      Line(points = {{-142, -70}, {-142, -75.5}, {-162, -75.5}, {-162, -67}}, color = {95, 95, 95}));
    connect(bodyCylinder5.frame_a, revolute3.frame_b) annotation (
      Line(points = {{-140, -35}, {-140, -29.5}, {-166, -29.5}, {-166, -55.25}, {-162, -55.25}, {-162, -67}}, color = {95, 95, 95}));
    connect(bodyCylinder8.frame_a, revolute4.frame_b) annotation (
      Line(points = {{-136, 158}, {-136, 157.5}, {-158, 157.5}, {-158, 127}}, color = {95, 95, 95}));
    connect(bodyCylinder7.frame_a, revolute4.frame_b) annotation (
      Line(points = {{-134, 114}, {-134, 110.375}, {-158, 110.375}, {-158, 127}}, color = {95, 95, 95}));
    connect(fixed.frame_b, fixedTranslation.frame_a) annotation (
      Line(points = {{-224, 46}, {-208, 46}}));
    connect(fixedTranslation.frame_b, revolute1.frame_a) annotation (
      Line(points = {{-188, 46}, {-188, 55}, {-178, 55}}));
    connect(fixedTranslation.frame_b, fixedTranslation1.frame_a) annotation (
      Line(points = {{-188, 46}, {-188, 32}}));
    connect(fixedTranslation1.frame_b, revolute2.frame_a) annotation (
      Line(points = {{-188, 12}, {-188, 11}, {-168, 11}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_b, revolute3.frame_a) annotation (
      Line(points = {{-190, -16}, {-191.25, -16}, {-191.25, -32}, {-192.5, -32}, {-192.5, -40}, {-184, -40}, {-184, -67}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_a, fixedTranslation1.frame_b) annotation (
      Line(points = {{-190, 4}, {-189, 4}, {-189, 12}, {-188, 12}}));
    connect(fixedTranslation3.frame_b, revolute5.frame_a) annotation (
      Line(points = {{-290, -126}, {-303, -126}, {-303, -125}, {-274, -125}}, color = {95, 95, 95}));
    connect(revolute5.frame_b, fixedTranslation4.frame_a) annotation (
      Line(points = {{-248, -125}, {-248, -137.5}, {-228, -137.5}, {-228, -124}}));
    connect(bodyCylinder10.frame_a, revolute6.frame_b) annotation (
      Line(points = {{-138, -148}, {-138, -151.813}, {-146, -151.813}, {-146, -149.625}, {-170, -149.625}, {-170, -121}}, color = {95, 95, 95}));
    connect(bodyCylinder9.frame_a, revolute6.frame_b) annotation (
      Line(points = {{-138, -110}, {-138, -106.75}, {-160, -106.75}, {-160, -105.5}, {-170, -105.5}, {-170, -121}}, color = {95, 95, 95}));
    connect(bodyCylinder3.frame_a, revolute1.frame_b) annotation (
      Line(points = {{-134, 86}, {-156, 86}, {-156, 55}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_b, fixedTranslation3.frame_a) annotation (
      Line(points = {{-190, -16}, {-190, -20}, {-310, -20}, {-310, -126}}, color = {95, 95, 95}));
    connect(bodyCylinder11.frame_a, revolute7.frame_b) annotation (
      Line(points={{-144,-226},{-144,-238.813},{-174,-238.813},{-174,-203}},          color = {95, 95, 95}));
    connect(bodyCylinder12.frame_a, revolute7.frame_b) annotation (
      Line(points = {{-146, -182}, {-146, -203.5}, {-174, -203.5}, {-174, -203}}, color = {95, 95, 95}));
    connect(fixedTranslation5.frame_b, revolute8.frame_a) annotation (
      Line(points = {{-300, -222}, {-300, -224}, {-286, -224}}, color = {95, 95, 95}));
    connect(fixedTranslation1.frame_b, fixedTranslation5.frame_a) annotation (
      Line(points = {{-188, 12}, {-320, 12}, {-320, -222}}, color = {95, 95, 95}));
    connect(fixed.frame_b, fixedTranslation7.frame_a) annotation (
      Line(points = {{-224, 46}, {-341.75, 46}, {-341.75, 176}, {-338.875, 176}, {-338.875, 212}, {-298, 212}}));
    connect(bodyCylinder13.frame_a, revolute9.frame_b) annotation (
      Line(points = {{-132, 194}, {-132, 185.188}, {-144, 185.188}, {-144, 184.375}, {-164, 184.375}, {-164, 217}}, color = {95, 95, 95}));
    connect(bodyCylinder14.frame_a, revolute9.frame_b) annotation (
      Line(points = {{-132, 230}, {-132, 240.5}, {-164, 240.5}, {-164, 217}}, color = {95, 95, 95}));
    connect(revolute9.frame_a, fixedTranslation9.frame_b) annotation (
      Line(points={{-186,217},{-198,217},{-198,216}},        color = {95, 95, 95}));
    connect(revolute10.frame_b, fixedTranslation9.frame_a) annotation (
      Line(points = {{-236, 218}, {-223, 218}, {-223, 216}, {-218, 216}}));
    connect(fixedTranslation7.frame_b, revolute10.frame_a) annotation (
      Line(points = {{-278, 212}, {-256, 212}, {-256, 218}}, color = {95, 95, 95}));
    connect(revolute8.frame_b, fixedTranslation6.frame_a) annotation (
      Line(points = {{-266, -224}, {-240.5, -224}, {-240.5, -202}, {-236, -202}}));
    connect(revolute7.frame_a, fixedTranslation6.frame_b) annotation (
      Line(points = {{-196, -203}, {-216, -203}, {-216, -202}}, color = {95, 95, 95}));
    connect(revolute6.frame_a, fixedTranslation4.frame_b) annotation (
      Line(points = {{-192, -121}, {-208, -121}, {-208, -124}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_b, fixedTranslation10.frame_a) annotation (
      Line(points = {{-190, -16}, {-190, -32}, {-198, -32}, {-198, -48}}, color = {95, 95, 95}));
    connect(fixedTranslation8.frame_a, revolute5.frame_b) annotation (
      Line(points={{-298,-84},{-298,-114},{-248,-114},{-248,-125}}));
    connect(springDamperParallel.frame_b, fixedTranslation8.frame_b) annotation (
      Line(points = {{-264, -68}, {-252, -68}, {-252, -84}, {-278, -84}}, color = {95, 95, 95}));
    connect(springDamperParallel.frame_a, fixedTranslation10.frame_b) annotation (
      Line(points = {{-264, -48}, {-218, -48}}, color = {95, 95, 95}));
    connect(fixedTranslation12.frame_a, fixedTranslation1.frame_b) annotation (
      Line(points = {{-364, -10}, {-220, -10}, {-220, 12}, {-188, 12}}));
    connect(fixedTranslation12.frame_b, springDamperParallel1.frame_a) annotation (
      Line(points = {{-384, -10}, {-401, -10}, {-401, -72}}, color = {95, 95, 95}));
    connect(fixedTranslation11.frame_a, revolute8.frame_b) annotation (
      Line(points = {{-378, -182}, {-266, -182}, {-266, -224}}, color = {95, 95, 95}));
    connect(fixedTranslation11.frame_b, springDamperParallel1.frame_b) annotation (
      Line(points = {{-398, -182}, {-401, -182}, {-401, -94}}, color = {95, 95, 95}));
    connect(fixedTranslation14.frame_a, fixedTranslation1.frame_b) annotation (
      Line(points = {{-340, 30}, {-188, 30}, {-188, 12}}, color = {95, 95, 95}));
    connect(fixedTranslation14.frame_b, springDamperParallel2.frame_a) annotation (
      Line(points = {{-360, 30}, {-410, 30}, {-410, 46}}, color = {95, 95, 95}));
    connect(revolute10.frame_b, fixedTranslation13.frame_a) annotation (
      Line(points = {{-236, 218}, {-224, 218}, {-224, 152}, {-250, 152}}));
    connect(fixedTranslation13.frame_b, springDamperParallel2.frame_b) annotation (
      Line(points = {{-250, 132}, {-410, 132}, {-410, 66}}));
    connect(bodyCylinder.frame_b, contact_sprocket.frame_b) annotation (
      Line(points={{-110,276},{16.2,276},{16.2,240.4}}));
    connect(contact_sprocket.frame_a, cinghia.frame_b) annotation (
      Line(points={{38,240.2},{136,240.2},{136,21.82}}, color = {95, 95, 95}, thickness = 0.5));
    connect(contact6.frame_a, cinghia.frame_b) annotation (
      Line(points={{-110,-321.36},{156,-321.36},{156,21.82},{136,21.82}},
                                                                       thickness = 0.5));
    connect(revolute7.frame_b, contact6.frame_b) annotation (
      Line(points={{-174,-203},{-228,-203},{-228,-320.72},{-179.76,-320.72}}, color = {95, 95, 95}));
    connect(contact5.frame_a, cinghia.frame_b) annotation (
      Line(points={{-18,-298.34},{136,-298.34},{136,21.82}},
                                                           color = {95, 95, 95}, thickness = 0.5));
    connect(revolute6.frame_b, contact5.frame_b) annotation (
      Line(points={{-170,-121},{-126,-121},{-126,-297.68},{-89.94,-297.68}}, color = {95, 95, 95}));
    connect(contact4.frame_a, cinghia.frame_b) annotation (
      Line(points={{-12,-209.48},{136,-209.48},{136,21.82}},
                                                           color = {95, 95, 95}, thickness = 0.5));
    connect(revolute3.frame_b, contact4.frame_b) annotation (
      Line(points={{-162,-67},{-92,-67},{-92,-208.96},{-68.68,-208.96}}, color = {95, 95, 95}));
    connect(contact3.frame_a, cinghia.frame_b) annotation (
      Line(points={{8,-134.46},{136,-134.46},{136,21.82}},
                                                         thickness = 0.5));
    connect(revolute2.frame_b, contact3.frame_b) annotation (
      Line(points={{-146,11},{-74,11},{-74,-133.92},{-50.86,-133.92}}));
    connect(contact2.frame_a, cinghia.frame_b) annotation (
      Line(points={{28,-43.32},{136,-43.32},{136,21.82}},
                                                        color = {95, 95, 95}, thickness = 0.5));
    connect(revolute1.frame_b, contact2.frame_b) annotation (
      Line(points={{-156,55},{-68,55},{-68,-42.64},{-46.12,-42.64}}, color = {95, 95, 95}));
    connect(contact.frame_a, cinghia.frame_b) annotation (
      Line(points={{40,174.6},{96,174.6},{96,21.82},{136,21.82}},thickness = 0.5));
    connect(revolute9.frame_b, contact.frame_b) annotation (
      Line(points={{-164,217},{-42,217},{-42,175.2},{-25.4,175.2}}));
    connect(contact1.frame_a, cinghia.frame_b) annotation (
      Line(points={{48,68.72},{72,68.72},{72,21.82},{136,21.82}},
                                                               color = {95, 95, 95}, thickness = 0.5));
    connect(revolute4.frame_b, contact1.frame_b) annotation (
      Line(points={{-158,127},{-44,127},{-44,69.44},{-30.48,69.44}}, color = {95, 95, 95}));
    connect(position.flange, revolute.axis) annotation (
      Line(points = {{-402, 300}, {-186, 300}, {-186, 280}}));
    connect(realExpression.y, position.phi_ref) annotation (
      Line(points={{-517,292},{-424,292},{-424,300}},        color = {0, 0, 127}));
    connect(prismatic.frame_b, revolute4.frame_a) annotation (Line(
        points={{-210,104},{-202,104},{-202,110},{-180,110},{-180,127}},
        color={95,95,95},
        thickness=0.5));
    connect(prismatic.frame_a, fixed.frame_b) annotation (Line(
        points={{-210,84},{-218,84},{-218,46},{-224,46}},
        color={95,95,95},
        thickness=0.5));
    connect(position1.flange, prismatic.axis) annotation (Line(points={{-240,100},
            {-228,100},{-228,102},{-216,102}}, color={0,127,0}));
    connect(realExpression1.y, position1.s_ref) annotation (Line(points={{-293,68},
            {-278,68},{-278,100},{-262,100}}, color={0,0,127}));
    annotation (
      Icon(coordinateSystem(grid = {2, 0})),
      experiment(StartTime = 0, StopTime = 20, Tolerance = 1e-03, Interval = 0.04));
  end track;

  model cinghia
    parameter Real sX=1.130312455;
    parameter Real sY=-0.2;
    parameter Real sZ=0.05;
    parameter Real l = 2.28;
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b[80] annotation (
      Placement(visible = true, transformation(origin = {-100, -3}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-100, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a[80] annotation (
      Placement(visible = true, transformation(origin = {100, -9}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {100, -8}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    AGILEX_BUNKER_PRO.element link2(
      startX=sX,
      startY=sY,
      startZ=sZ,
      thZ=Modelica.Constants.pi/2,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link21(
      startX=link2.startX + link2.bodyShape.r[1]*cos(link2.thZ),
      startY=link2.startY + link2.bodyShape.r[1]*sin(link2.thZ),
      startZ=sZ,
      thZ=link2.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link22(
      startX=link21.startX + link21.bodyShape.r[1]*cos(link21.thZ),
      startY=link21.startY + link21.bodyShape.r[1]*sin(link21.thZ),
      startZ=sZ,
      thZ=link21.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link23(
      startX=link22.startX + link22.bodyShape.r[1]*cos(link22.thZ),
      startY=link22.startY + link22.bodyShape.r[1]*sin(link22.thZ),
      startZ=sZ,
      thZ=link22.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link24(
      startX=link23.startX + link23.bodyShape.r[1]*cos(link23.thZ),
      startY=link23.startY + link23.bodyShape.r[1]*sin(link23.thZ),
      startZ=sZ,
      thZ=link23.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link25(
      startX=link24.startX + link24.bodyShape.r[1]*cos(link24.thZ),
      startY=link24.startY + link24.bodyShape.r[1]*sin(link24.thZ),
      startZ=sZ,
      thZ=link24.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link26(
      startX=link25.startX + link25.bodyShape.r[1]*cos(link25.thZ),
      startY=link25.startY + link25.bodyShape.r[1]*sin(link25.thZ),
      startZ=sZ,
      thZ=link25.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link27(
      startX=link26.startX + link26.bodyShape.r[1]*cos(link26.thZ),
      startY=link26.startY + link26.bodyShape.r[1]*sin(link26.thZ),
      startZ=sZ,
      thZ=link26.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link28(
      startX=link27.startX + link27.bodyShape.r[1]*cos(link27.thZ),
      startY=link27.startY + link27.bodyShape.r[1]*sin(link27.thZ),
      startZ=sZ,
      thZ=link27.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link29(
      startX=link28.startX + link28.bodyShape.r[1]*cos(link28.thZ),
      startY=link28.startY + link28.bodyShape.r[1]*sin(link28.thZ),
      startZ=sZ,
      thZ=link28.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link210(
      startX=link29.startX + link29.bodyShape.r[1]*cos(link29.thZ),
      startY=link29.startY + link29.bodyShape.r[1]*sin(link29.thZ),
      startZ=sZ,
      thZ=link29.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link211(
      startX=link210.startX + link210.bodyShape.r[1]*cos(link210.thZ),
      startY=link210.startY + link210.bodyShape.r[1]*sin(link210.thZ),
      startZ=sZ,
      thZ=link210.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link212(
      startX=link211.startX + link211.bodyShape.r[1]*cos(link211.thZ),
      startY=link211.startY + link211.bodyShape.r[1]*sin(link211.thZ),
      startZ=sZ,
      thZ=link211.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link213(
      startX=link212.startX + link212.bodyShape.r[1]*cos(link212.thZ),
      startY=link212.startY + link212.bodyShape.r[1]*sin(link212.thZ),
      startZ=sZ,
      thZ=link212.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link214(
      startX=link213.startX + link213.bodyShape.r[1]*cos(link213.thZ),
      startY=link213.startY + link213.bodyShape.r[1]*sin(link213.thZ),
      startZ=sZ,
      thZ=link213.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link215(
      startX=link214.startX + link214.bodyShape.r[1]*cos(link214.thZ),
      startY=link214.startY + link214.bodyShape.r[1]*sin(link214.thZ),
      startZ=sZ,
      thZ=link214.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link216(
      startX=link215.startX + link215.bodyShape.r[1]*cos(link215.thZ),
      startY=link215.startY + link215.bodyShape.r[1]*sin(link215.thZ),
      startZ=sZ,
      thZ=link215.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link217(
      startX=link216.startX + link216.bodyShape.r[1]*cos(link216.thZ),
      startY=link216.startY + link216.bodyShape.r[1]*sin(link216.thZ),
      startZ=sZ,
      thZ=link216.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link218(
      startX=link217.startX + link217.bodyShape.r[1]*cos(link217.thZ),
      startY=link217.startY + link217.bodyShape.r[1]*sin(link217.thZ),
      startZ=sZ,
      thZ=link217.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link219(
      startX=link218.startX + link218.bodyShape.r[1]*cos(link218.thZ),
      startY=link218.startY + link218.bodyShape.r[1]*sin(link218.thZ),
      startZ=sZ,
      thZ=link218.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    Modelica.Mechanics.MultiBody.Joints.Prismatic prismatic(animation=true,    useAxisFlange = true) annotation (
      Placement(visible = true, transformation(origin={-4,-54},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Translational.Sources.Position position(exact=false)
                                                             annotation (
      Placement(visible = true, transformation(origin = {-12, -15}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp(
      duration=3,
      height=-l*2/Modelica.Constants.pi,
      offset=l*2/Modelica.Constants.pi,
     startTime=0)                                                                                                                                       annotation (
     Placement(visible = true, transformation(origin={-70,-30},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    AGILEX_BUNKER_PRO.element link220(
      startX=link219.startX + link219.bodyShape.r[1]*cos(link219.thZ),
      startY=link219.startY + link219.bodyShape.r[1]*sin(link219.thZ),
      startZ=sZ,
      thZ=link219.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link221(
      startX=link220.startX + link220.bodyShape.r[1]*cos(link220.thZ),
      startY=link220.startY + link220.bodyShape.r[1]*sin(link220.thZ),
      startZ=sZ,
      thZ=link220.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link222(
      startX=link221.startX + link221.bodyShape.r[1]*cos(link221.thZ),
      startY=link221.startY + link221.bodyShape.r[1]*sin(link221.thZ),
      startZ=sZ,
      thZ=link221.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link223(
      startX=link222.startX + link222.bodyShape.r[1]*cos(link222.thZ),
      startY=link222.startY + link222.bodyShape.r[1]*sin(link222.thZ),
      startZ=sZ,
      thZ=link222.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link224(
      startX=link223.startX + link223.bodyShape.r[1]*cos(link223.thZ),
      startY=link223.startY + link223.bodyShape.r[1]*sin(link223.thZ),
      startZ=sZ,
      thZ=link223.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link225(
      startX=link224.startX + link224.bodyShape.r[1]*cos(link224.thZ),
      startY=link224.startY + link224.bodyShape.r[1]*sin(link224.thZ),
      startZ=sZ,
      thZ=link224.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
      AGILEX_BUNKER_PRO.element link226(
      startX=link225.startX + link225.bodyShape.r[1]*cos(link225.thZ),
      startY=link225.startY + link225.bodyShape.r[1]*sin(link225.thZ),
      startZ=sZ,
      thZ=link225.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link227(
      startX=link226.startX + link226.bodyShape.r[1]*cos(link226.thZ),
      startY=link226.startY + link226.bodyShape.r[1]*sin(link226.thZ),
      startZ=sZ,
      thZ=link226.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link228(
      startX=link227.startX + link227.bodyShape.r[1]*cos(link227.thZ),
      startY=link227.startY + link227.bodyShape.r[1]*sin(link227.thZ),
      startZ=sZ,
      thZ=link227.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link229(
      startX=link228.startX + link228.bodyShape.r[1]*cos(link228.thZ),
      startY=link228.startY + link228.bodyShape.r[1]*sin(link228.thZ),
      startZ=sZ,
      thZ=link228.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
       AGILEX_BUNKER_PRO.element link230(
      startX=link229.startX + link229.bodyShape.r[1]*cos(link229.thZ),
      startY=link229.startY + link229.bodyShape.r[1]*sin(link229.thZ),
      startZ=sZ,
      thZ=link229.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
      AGILEX_BUNKER_PRO.element link231(
      startX=link230.startX + link230.bodyShape.r[1]*cos(link230.thZ),
      startY=link230.startY + link230.bodyShape.r[1]*sin(link230.thZ),
      startZ=sZ,
      thZ=link230.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
      AGILEX_BUNKER_PRO.element link232(
      startX=link231.startX + link231.bodyShape.r[1]*cos(link231.thZ),
      startY=link231.startY + link231.bodyShape.r[1]*sin(link231.thZ),
      startZ=sZ,
      thZ=link231.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
      AGILEX_BUNKER_PRO.element link233(
      startX=link232.startX + link232.bodyShape.r[1]*cos(link232.thZ),
      startY=link232.startY + link232.bodyShape.r[1]*sin(link232.thZ),
      startZ=sZ,
      thZ=link232.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
       AGILEX_BUNKER_PRO.element link234(
      startX=link233.startX + link233.bodyShape.r[1]*cos(link233.thZ),
      startY=link233.startY + link233.bodyShape.r[1]*sin(link233.thZ),
      startZ=sZ,
      thZ=link233.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
      AGILEX_BUNKER_PRO.element link235(
      startX=link234.startX + link234.bodyShape.r[1]*cos(link234.thZ),
      startY=link234.startY + link234.bodyShape.r[1]*sin(link234.thZ),
      startZ=sZ,
      thZ=link234.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link236(
      startX=link235.startX + link235.bodyShape.r[1]*cos(link235.thZ),
      startY=link235.startY + link235.bodyShape.r[1]*sin(link235.thZ),
      startZ=sZ,
      thZ=link235.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
      AGILEX_BUNKER_PRO.element link237(
      startX=link236.startX + link236.bodyShape.r[1]*cos(link236.thZ),
      startY=link236.startY + link236.bodyShape.r[1]*sin(link236.thZ),
      startZ=sZ,
      thZ=link236.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link238(
      startX=link237.startX + link237.bodyShape.r[1]*cos(link237.thZ),
      startY=link237.startY + link237.bodyShape.r[1]*sin(link237.thZ),
      startZ=sZ,
      thZ=link237.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link239(
      startX=link238.startX + link238.bodyShape.r[1]*cos(link238.thZ),
      startY=link238.startY + link238.bodyShape.r[1]*sin(link238.thZ),
      startZ=sZ,
      thZ=link238.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link240(
      startX=link239.startX + link239.bodyShape.r[1]*cos(link239.thZ),
      startY=link239.startY + link239.bodyShape.r[1]*sin(link239.thZ),
      startZ=sZ,
      thZ=link239.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link241(
      startX=link240.startX + link240.bodyShape.r[1]*cos(link240.thZ),
      startY=link240.startY + link240.bodyShape.r[1]*sin(link240.thZ),
      startZ=sZ,
      thZ=link240.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link242(
      startX=link241.startX + link241.bodyShape.r[1]*cos(link241.thZ),
      startY=link241.startY + link241.bodyShape.r[1]*sin(link241.thZ),
      startZ=sZ,
      thZ=link241.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link243(
      startX=link242.startX + link242.bodyShape.r[1]*cos(link242.thZ),
      startY=link242.startY + link242.bodyShape.r[1]*sin(link242.thZ),
      startZ=sZ,
      thZ=link242.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link244(
      startX=link243.startX + link243.bodyShape.r[1]*cos(link243.thZ),
      startY=link243.startY + link243.bodyShape.r[1]*sin(link243.thZ),
      startZ=sZ,
      thZ=link243.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link245(
      startX=link244.startX + link244.bodyShape.r[1]*cos(link244.thZ),
      startY=link244.startY + link244.bodyShape.r[1]*sin(link244.thZ),
      startZ=sZ,
      thZ=link244.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link246(
      startX=link245.startX + link245.bodyShape.r[1]*cos(link245.thZ),
      startY=link245.startY + link245.bodyShape.r[1]*sin(link245.thZ),
      startZ=sZ,
      thZ=link245.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link247(
      startX=link246.startX + link246.bodyShape.r[1]*cos(link246.thZ),
      startY=link246.startY + link246.bodyShape.r[1]*sin(link246.thZ),
      startZ=sZ,
      thZ=link246.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link248(
      startX=link247.startX + link247.bodyShape.r[1]*cos(link247.thZ),
      startY=link247.startY + link247.bodyShape.r[1]*sin(link247.thZ),
      startZ=sZ,
      thZ=link247.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
      AGILEX_BUNKER_PRO.element link249(
      startX=link248.startX + link248.bodyShape.r[1]*cos(link248.thZ),
      startY=link248.startY + link248.bodyShape.r[1]*sin(link248.thZ),
      startZ=sZ,
      thZ=link248.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link250(
      startX=link249.startX + link249.bodyShape.r[1]*cos(link249.thZ),
      startY=link249.startY + link249.bodyShape.r[1]*sin(link249.thZ),
      startZ=sZ,
      thZ=link249.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link251(
      startX=link250.startX + link250.bodyShape.r[1]*cos(link250.thZ),
      startY=link250.startY + link250.bodyShape.r[1]*sin(link250.thZ),
      startZ=sZ,
      thZ=link250.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link252(
      startX=link251.startX + link251.bodyShape.r[1]*cos(link251.thZ),
      startY=link251.startY + link251.bodyShape.r[1]*sin(link251.thZ),
      startZ=sZ,
      thZ=link251.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link253(
      startX=link252.startX + link252.bodyShape.r[1]*cos(link252.thZ),
      startY=link252.startY + link252.bodyShape.r[1]*sin(link252.thZ),
      startZ=sZ,
      thZ=link252.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
      AGILEX_BUNKER_PRO.element link254(
      startX=link253.startX + link253.bodyShape.r[1]*cos(link253.thZ),
      startY=link253.startY + link253.bodyShape.r[1]*sin(link253.thZ),
      startZ=sZ,
      thZ=link253.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
      AGILEX_BUNKER_PRO.element link255(
      startX=link254.startX + link254.bodyShape.r[1]*cos(link254.thZ),
      startY=link254.startY + link254.bodyShape.r[1]*sin(link254.thZ),
      startZ=sZ,
      thZ=link254.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
      AGILEX_BUNKER_PRO.element link256(
      startX=link255.startX + link255.bodyShape.r[1]*cos(link255.thZ),
      startY=link255.startY + link255.bodyShape.r[1]*sin(link255.thZ),
      startZ=sZ,
      thZ=link255.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link257(
      startX=link256.startX + link256.bodyShape.r[1]*cos(link256.thZ),
      startY=link256.startY + link256.bodyShape.r[1]*sin(link256.thZ),
      startZ=sZ,
      thZ=link256.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link258(
      startX=link257.startX + link257.bodyShape.r[1]*cos(link257.thZ),
      startY=link257.startY + link257.bodyShape.r[1]*sin(link257.thZ),
      startZ=sZ,
      thZ=link257.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link259(
      startX=link258.startX + link258.bodyShape.r[1]*cos(link258.thZ),
      startY=link258.startY + link258.bodyShape.r[1]*sin(link258.thZ),
      startZ=sZ,
      thZ=link258.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link260(
      startX=link259.startX + link259.bodyShape.r[1]*cos(link259.thZ),
      startY=link259.startY + link259.bodyShape.r[1]*sin(link259.thZ),
      startZ=sZ,
      thZ=link259.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
      AGILEX_BUNKER_PRO.element link261(
      startX=link260.startX + link260.bodyShape.r[1]*cos(link260.thZ),
      startY=link260.startY + link260.bodyShape.r[1]*sin(link260.thZ),
      startZ=sZ,
      thZ=link260.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link262(
      startX=link261.startX + link261.bodyShape.r[1]*cos(link261.thZ),
      startY=link261.startY + link261.bodyShape.r[1]*sin(link261.thZ),
      startZ=sZ,
      thZ=link261.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link263(
      startX=link262.startX + link262.bodyShape.r[1]*cos(link262.thZ),
      startY=link262.startY + link262.bodyShape.r[1]*sin(link262.thZ),
      startZ=sZ,
      thZ=link262.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link264(
      startX=link263.startX + link263.bodyShape.r[1]*cos(link263.thZ),
      startY=link263.startY + link263.bodyShape.r[1]*sin(link263.thZ),
      startZ=sZ,
      thZ=link263.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link265(
      startX=link264.startX + link264.bodyShape.r[1]*cos(link264.thZ),
      startY=link264.startY + link264.bodyShape.r[1]*sin(link264.thZ),
      startZ=sZ,
      thZ=link264.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link266(
      startX=link265.startX + link265.bodyShape.r[1]*cos(link265.thZ),
      startY=link265.startY + link265.bodyShape.r[1]*sin(link265.thZ),
      startZ=sZ,
      thZ=link265.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link267(
      startX=link266.startX + link266.bodyShape.r[1]*cos(link266.thZ),
      startY=link266.startY + link266.bodyShape.r[1]*sin(link266.thZ),
      startZ=sZ,
      thZ=link266.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link268(
      startX=link267.startX + link267.bodyShape.r[1]*cos(link267.thZ),
      startY=link267.startY + link267.bodyShape.r[1]*sin(link267.thZ),
      startZ=sZ,
      thZ=link267.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link269(
      startX=link268.startX + link268.bodyShape.r[1]*cos(link268.thZ),
      startY=link268.startY + link268.bodyShape.r[1]*sin(link268.thZ),
      startZ=sZ,
      thZ=link268.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    AGILEX_BUNKER_PRO.element link270(
      startX=link269.startX + link269.bodyShape.r[1]*cos(link269.thZ),
      startY=link269.startY + link269.bodyShape.r[1]*sin(link269.thZ),
      startZ=sZ,
      thZ=link269.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link271(
      startX=link270.startX + link270.bodyShape.r[1]*cos(link270.thZ),
      startY=link270.startY + link270.bodyShape.r[1]*sin(link270.thZ),
      startZ=sZ,
      thZ=link270.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link272(
      startX=link271.startX + link271.bodyShape.r[1]*cos(link271.thZ),
      startY=link271.startY + link271.bodyShape.r[1]*sin(link271.thZ),
      startZ=sZ,
      thZ=link271.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link273(
      startX=link272.startX + link272.bodyShape.r[1]*cos(link272.thZ),
      startY=link272.startY + link272.bodyShape.r[1]*sin(link272.thZ),
      startZ=sZ,
      thZ=link272.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link274(
      startX=link273.startX + link273.bodyShape.r[1]*cos(link273.thZ),
      startY=link273.startY + link273.bodyShape.r[1]*sin(link273.thZ),
      startZ=sZ,
      thZ=link273.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link275(
      startX=link274.startX + link274.bodyShape.r[1]*cos(link274.thZ),
      startY=link274.startY + link274.bodyShape.r[1]*sin(link274.thZ),
      startZ=sZ,
      thZ=link274.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link276(
      startX=link275.startX + link275.bodyShape.r[1]*cos(link275.thZ),
      startY=link275.startY + link275.bodyShape.r[1]*sin(link275.thZ),
      startZ=sZ,
      thZ=link275.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
      AGILEX_BUNKER_PRO.element link277(
      startX=link276.startX + link276.bodyShape.r[1]*cos(link276.thZ),
      startY=link276.startY + link276.bodyShape.r[1]*sin(link276.thZ),
      startZ=sZ,
      thZ=link276.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
      AGILEX_BUNKER_PRO.element link278(
      startX=link277.startX + link277.bodyShape.r[1]*cos(link277.thZ),
      startY=link277.startY + link277.bodyShape.r[1]*sin(link277.thZ),
      startZ=sZ,
      thZ=link277.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
     AGILEX_BUNKER_PRO.element link279(
      startX=link278.startX + link278.bodyShape.r[1]*cos(link278.thZ),
      startY=link278.startY + link278.bodyShape.r[1]*sin(link278.thZ),
      startZ=sZ,
      thZ=link278.thZ + Modelica.Constants.pi/80,
      l=l) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
    Modelica.Mechanics.MultiBody.Parts.FixedRotation fixedRotation(n(
          displayUnit="1") = {0,0,1},                              angle=90)
      annotation (Placement(transformation(extent={{28,-66},{48,-48}})));
    Modelica.Mechanics.MultiBody.Parts.FixedRotation fixedRotation1(n(
          displayUnit="1") = {0,0,1},                               angle=90)
      annotation (Placement(transformation(extent={{-98,-66},{-78,-48}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute
      annotation (Placement(transformation(extent={{-56,-66},{-36,-45}})));
  equation
    connect(link2.frame_b, link21.frame_a);
    connect(link21.frame_b, link22.frame_a);
    connect(link22.frame_b, link23.frame_a);
    connect(link23.frame_b, link24.frame_a);
    connect(link24.frame_b, link25.frame_a);
    connect(link25.frame_b, link26.frame_a);
    connect(link26.frame_b, link27.frame_a);
    connect(link27.frame_b, link28.frame_a);
    connect(link28.frame_b, link29.frame_a);
    connect(link29.frame_b, link210.frame_a);
    connect(link210.frame_b, link211.frame_a);
    connect(link211.frame_b, link212.frame_a);
    connect(link212.frame_b, link213.frame_a);
    connect(link213.frame_b, link214.frame_a);
    connect(link214.frame_b, link215.frame_a);
    connect(link215.frame_b, link216.frame_a);
    connect(link216.frame_b, link217.frame_a);
    connect(link217.frame_b, link218.frame_a);
    connect(link218.frame_b, link219.frame_a);
    connect(link219.frame_b, link220.frame_a);
    connect(link220.frame_b, link221.frame_a);
    connect(link221.frame_b, link222.frame_a);
    connect(link222.frame_b, link223.frame_a);
    connect(link223.frame_b, link224.frame_a);
    connect(link224.frame_b, link225.frame_a);
    connect(link225.frame_b, link226.frame_a);
    connect(link226.frame_b, link227.frame_a);
    connect(link227.frame_b, link228.frame_a);
    connect(link228.frame_b, link229.frame_a);
    connect(link229.frame_b, link230.frame_a);
    connect(link230.frame_b, link231.frame_a);
    connect(link231.frame_b, link232.frame_a);
    connect(link232.frame_b, link233.frame_a);
    connect(link233.frame_b, link234.frame_a);
    connect(link234.frame_b, link235.frame_a);
    connect(link235.frame_b, link236.frame_a);
    connect(link236.frame_b, link237.frame_a);
    connect(link237.frame_b, link238.frame_a);
    connect(link238.frame_b, link239.frame_a);
    connect(link239.frame_b, link240.frame_a);
    connect(link240.frame_b, link241.frame_a);
    connect(link241.frame_b, link242.frame_a);
    connect(link242.frame_b, link243.frame_a);
    connect(link243.frame_b, link244.frame_a);
    connect(link244.frame_b, link245.frame_a);
    connect(link245.frame_b, link246.frame_a);
    connect(link246.frame_b, link247.frame_a);
    connect(link247.frame_b, link248.frame_a);
    connect(link248.frame_b, link249.frame_a);
    connect(link249.frame_b, link250.frame_a);
    connect(link250.frame_b, link251.frame_a);
    connect(link251.frame_b, link252.frame_a);
    connect(link252.frame_b, link253.frame_a);
    connect(link253.frame_b, link254.frame_a);
    connect(link254.frame_b, link255.frame_a);
    connect(link255.frame_b, link256.frame_a);
    connect(link256.frame_b, link257.frame_a);
    connect(link257.frame_b, link258.frame_a);
    connect(link258.frame_b, link259.frame_a);
    connect(link259.frame_b, link260.frame_a);
    connect(link260.frame_b, link261.frame_a);
    connect(link261.frame_b, link262.frame_a);
    connect(link262.frame_b, link263.frame_a);
    connect(link263.frame_b, link264.frame_a);
    connect(link264.frame_b, link265.frame_a);
    connect(link265.frame_b, link266.frame_a);
    connect(link266.frame_b, link267.frame_a);
    connect(link267.frame_b, link268.frame_a);
    connect(link268.frame_b, link269.frame_a);
    connect(link269.frame_b, link270.frame_a);
    connect(link270.frame_b, link271.frame_a);
    connect(link271.frame_b, link272.frame_a);
    connect(link272.frame_b, link273.frame_a);
    connect(link273.frame_b, link274.frame_a);
    connect(link274.frame_b, link275.frame_a);
    connect(link275.frame_b, link276.frame_a);
    connect(link276.frame_b, link277.frame_a);
    connect(link277.frame_b, link278.frame_a);
    connect(link278.frame_b, link279.frame_a);
    connect(link279.frame_b, fixedRotation1.frame_a);
    connect(fixedRotation.frame_b, link2.frame_a);
    connect(link2.frame_b1, frame_b[1]);
    connect(link21.frame_b1, frame_b[2]);
    connect(link22.frame_b1, frame_b[3]);
    connect(link23.frame_b1, frame_b[4]);
    connect(link24.frame_b1, frame_b[5]);
    connect(link25.frame_b1, frame_b[6]);
    connect(link26.frame_b1, frame_b[7]);
    connect(link27.frame_b1, frame_b[8]);
    connect(link28.frame_b1, frame_b[9]);
    connect(link29.frame_b1, frame_b[10]);
    connect(link210.frame_b1, frame_b[11]);
    connect(link211.frame_b1, frame_b[12]);
    connect(link212.frame_b1, frame_b[13]);
    connect(link213.frame_b1, frame_b[14]);
    connect(link214.frame_b1, frame_b[15]);
    connect(link215.frame_b1, frame_b[16]);
    connect(link216.frame_b1, frame_b[17]);
    connect(link217.frame_b1, frame_b[18]);
    connect(link218.frame_b1, frame_b[19]);
    connect(link219.frame_b1, frame_b[20]);
    connect(link220.frame_b1, frame_b[21]);
    connect(link221.frame_b1, frame_b[22]);
    connect(link222.frame_b1, frame_b[23]);
    connect(link223.frame_b1, frame_b[24]);
    connect(link224.frame_b1, frame_b[25]);
    connect(link225.frame_b1, frame_b[26]);
    connect(link226.frame_b1, frame_b[27]);
    connect(link227.frame_b1, frame_b[28]);
    connect(link228.frame_b1, frame_b[29]);
    connect(link229.frame_b1, frame_b[30]);
    connect(link230.frame_b1, frame_b[31]);
    connect(link231.frame_b1, frame_b[32]);
    connect(link232.frame_b1, frame_b[33]);
    connect(link233.frame_b1, frame_b[34]);
    connect(link234.frame_b1, frame_b[35]);
    connect(link235.frame_b1, frame_b[36]);
    connect(link236.frame_b1, frame_b[37]);
    connect(link237.frame_b1, frame_b[38]);
    connect(link238.frame_b1, frame_b[39]);
    connect(link239.frame_b1, frame_b[40]);
    connect(link240.frame_b1, frame_b[41]);
    connect(link241.frame_b1, frame_b[42]);
    connect(link242.frame_b1, frame_b[43]);
    connect(link243.frame_b1, frame_b[44]);
    connect(link244.frame_b1, frame_b[45]);
    connect(link245.frame_b1, frame_b[46]);
    connect(link246.frame_b1, frame_b[47]);
    connect(link247.frame_b1, frame_b[48]);
    connect(link248.frame_b1, frame_b[49]);
    connect(link249.frame_b1, frame_b[50]);
    connect(link250.frame_b1, frame_b[51]);
    connect(link251.frame_b1, frame_b[52]);
    connect(link252.frame_b1, frame_b[53]);
    connect(link253.frame_b1, frame_b[54]);
    connect(link254.frame_b1, frame_b[55]);
    connect(link255.frame_b1, frame_b[56]);
    connect(link256.frame_b1, frame_b[57]);
    connect(link257.frame_b1, frame_b[58]);
    connect(link258.frame_b1, frame_b[59]);
    connect(link259.frame_b1, frame_b[60]);
    connect(link260.frame_b1, frame_b[61]);
    connect(link261.frame_b1, frame_b[62]);
    connect(link262.frame_b1, frame_b[63]);
    connect(link263.frame_b1, frame_b[64]);
    connect(link264.frame_b1, frame_b[65]);
    connect(link265.frame_b1, frame_b[66]);
    connect(link266.frame_b1, frame_b[67]);
    connect(link267.frame_b1, frame_b[68]);
    connect(link268.frame_b1, frame_b[69]);
    connect(link269.frame_b1, frame_b[70]);
    connect(link270.frame_b1, frame_b[71]);
    connect(link271.frame_b1, frame_b[72]);
    connect(link272.frame_b1, frame_b[73]);
    connect(link273.frame_b1, frame_b[74]);
    connect(link274.frame_b1, frame_b[75]);
    connect(link275.frame_b1, frame_b[76]);
    connect(link276.frame_b1, frame_b[77]);
    connect(link277.frame_b1, frame_b[78]);
    connect(link278.frame_b1, frame_b[79]);
    connect(link279.frame_b1, frame_b[80]);

    connect(link2.frame_a1, frame_a[1]);
    connect(link21.frame_a1, frame_a[2]);
    connect(link22.frame_a1, frame_a[3]);
    connect(link23.frame_a1, frame_a[4]);
    connect(link24.frame_a1, frame_a[5]);
    connect(link25.frame_a1, frame_a[6]);
    connect(link26.frame_a1, frame_a[7]);
    connect(link27.frame_a1, frame_a[8]);
    connect(link28.frame_a1, frame_a[9]);
    connect(link29.frame_a1, frame_a[10]);
    connect(link210.frame_a1, frame_a[11]);
    connect(link211.frame_a1, frame_a[12]);
    connect(link212.frame_a1, frame_a[13]);
    connect(link213.frame_a1, frame_a[14]);
    connect(link214.frame_a1, frame_a[15]);
    connect(link215.frame_a1, frame_a[16]);
    connect(link216.frame_a1, frame_a[17]);
    connect(link217.frame_a1, frame_a[18]);
    connect(link218.frame_a1, frame_a[19]);
    connect(link219.frame_a1, frame_a[20]);
    connect(link220.frame_a1, frame_a[21]);
    connect(link221.frame_a1, frame_a[22]);
    connect(link222.frame_a1, frame_a[23]);
    connect(link223.frame_a1, frame_a[24]);
    connect(link224.frame_a1, frame_a[25]);
    connect(link225.frame_a1, frame_a[26]);
    connect(link226.frame_a1, frame_a[27]);
    connect(link227.frame_a1, frame_a[28]);
    connect(link228.frame_a1, frame_a[29]);
    connect(link229.frame_a1, frame_a[30]);
    connect(link230.frame_a1, frame_a[31]);
    connect(link231.frame_a1, frame_a[32]);
    connect(link232.frame_a1, frame_a[33]);
    connect(link233.frame_a1, frame_a[34]);
    connect(link234.frame_a1, frame_a[35]);
    connect(link235.frame_a1, frame_a[36]);
    connect(link236.frame_a1, frame_a[37]);
    connect(link237.frame_a1, frame_a[38]);
    connect(link238.frame_a1, frame_a[39]);
    connect(link239.frame_a1, frame_a[40]);
    connect(link240.frame_a1, frame_a[41]);
    connect(link241.frame_a1, frame_a[42]);
    connect(link242.frame_a1, frame_a[43]);
    connect(link243.frame_a1, frame_a[44]);
    connect(link244.frame_a1, frame_a[45]);
    connect(link245.frame_a1, frame_a[46]);
    connect(link246.frame_a1, frame_a[47]);
    connect(link247.frame_a1, frame_a[48]);
    connect(link248.frame_a1, frame_a[49]);
    connect(link249.frame_a1, frame_a[50]);
    connect(link250.frame_a1, frame_a[51]);
    connect(link251.frame_a1, frame_a[52]);
    connect(link252.frame_a1, frame_a[53]);
    connect(link253.frame_a1, frame_a[54]);
    connect(link254.frame_a1, frame_a[55]);
    connect(link255.frame_a1, frame_a[56]);
    connect(link256.frame_a1, frame_a[57]);
    connect(link257.frame_a1, frame_a[58]);
    connect(link258.frame_a1, frame_a[59]);
    connect(link259.frame_a1, frame_a[60]);
    connect(link260.frame_a1, frame_a[61]);
    connect(link261.frame_a1, frame_a[62]);
    connect(link262.frame_a1, frame_a[63]);
    connect(link263.frame_a1, frame_a[64]);
    connect(link264.frame_a1, frame_a[65]);
    connect(link265.frame_a1, frame_a[66]);
    connect(link266.frame_a1, frame_a[67]);
    connect(link267.frame_a1, frame_a[68]);
    connect(link268.frame_a1, frame_a[69]);
    connect(link269.frame_a1, frame_a[70]);
    connect(link270.frame_a1, frame_a[71]);
    connect(link271.frame_a1, frame_a[72]);
    connect(link272.frame_a1, frame_a[73]);
    connect(link273.frame_a1, frame_a[74]);
    connect(link274.frame_a1, frame_a[75]);
    connect(link275.frame_a1, frame_a[76]);
    connect(link276.frame_a1, frame_a[77]);
    connect(link277.frame_a1, frame_a[78]);
    connect(link278.frame_a1, frame_a[79]);
    connect(link279.frame_a1, frame_a[80]);
    connect(ramp.y, position.s_ref) annotation (Line(points={{-59,-30},{-32,-30},
            {-32,-15},{-24,-15}}, color={0,0,127}));
    connect(position.flange, prismatic.axis) annotation (Line(points={{-2,-15},{16,
            -15},{16,-48},{4,-48}},     color={0,127,0}));
    connect(prismatic.frame_b, fixedRotation.frame_a) annotation (Line(
        points={{6,-54},{6,-57},{28,-57}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedRotation1.frame_b, revolute.frame_a) annotation (Line(
        points={{-78,-57},{-56,-57},{-56,-55.5}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute.frame_b, prismatic.frame_a) annotation (Line(
        points={{-36,-55.5},{-22,-55.5},{-22,-54},{-14,-54}},
        color={95,95,95},
        thickness=0.5));
    annotation (
      Icon(coordinateSystem(grid = {2, 0})),
      Diagram(coordinateSystem(grid = {2, 3})));
  end cinghia;

  model contact_sprocket
    parameter Modelica.Units.SI.Length R = 0.08;
    AGILEX_BUNKER_PRO.contact_force_new contact_force[80](R=fill(R, 80))
      annotation (Placement(visible=true, transformation(
          origin={-4,3},
          extent={{10,-10},{-10,10}},
          rotation=0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b annotation (
      Placement(visible = true, transformation(origin = {-118, 3}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-118, 4}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a[80] annotation (
      Placement(visible = true, transformation(origin = {100, 3}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {100, 2}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
  equation
    for n in 1:80 loop
      connect(contact_force[n].frame_b, frame_b);
      connect(contact_force[n].frame_a, frame_a[n]);
    end for;
    annotation (
      Icon(coordinateSystem(grid = {2, 0})),
      Diagram(coordinateSystem(grid = {2, 3})));
  end contact_sprocket;

  model contact
    parameter Modelica.Units.SI.Length R = 0.08;
    parameter Real lar = 0.15;
    parameter Real L = 0.048;
    //parameter Real Drx = 10 ^ 2;
    //parameter Real Kry = 7*10 ^ 5;
    //parameter Real Dry = 10 ^ 2;
    AGILEX_BUNKER_PRO.contact_force_idler_new contact_force[80](
      R=fill(R, 80),
      each lar=lar,
      each L=L) annotation (Placement(visible=true, transformation(
          origin={-4,3},
          extent={{10,-10},{-10,10}},
          rotation=0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b annotation (
      Placement(visible = true, transformation(origin = {-118, 3}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-118, 4}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a[80] annotation (
      Placement(visible = true, transformation(origin = {100, 3}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {100, 2}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
  equation
    for n in 1:80 loop
      connect(contact_force[n].frame_b, frame_b);
      connect(contact_force[n].frame_a, frame_a[n]);
    end for;
    annotation (
      Icon(coordinateSystem(grid = {2, 0})),
      Diagram(coordinateSystem(grid = {2, 3})));
  end contact;

  model trackfinal
    extends Modelica.Icons.Example;
    inner Modelica.Mechanics.MultiBody.World world annotation (
      Placement(visible = true, transformation(origin = {-368, 186}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder(diameter = 0.16, length = 0.1, lengthDirection = {0, 0, 1}, r = {0, 0, 0.05}) annotation (
      Placement(visible = true, transformation(origin = {-120, 276}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute(animation = false, useAxisFlange = true) annotation (
      Placement(visible = true, transformation(origin = {-186, 270}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder4(diameter = 0.088, length = 0.075, r = {0, 0, 0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-124, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute1 annotation (
      Placement(visible = true, transformation(origin = {-167, 55}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder1(diameter = 0.088, length = 0.075, r = {0, 0, -0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-124, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder2(diameter = 0.088, length = 0.075, r = {0, 0, 0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-126, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute2 annotation (
      Placement(visible = true, transformation(origin = {-157, 11}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder5(diameter = 0.088, length = 0.075, r = {0, 0, -0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-129, -35}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder6(diameter = 0.088, length = 0.075, r = {0, 0, 0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-132, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute3 annotation (
      Placement(visible = true, transformation(origin = {-173, -67}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder7(diameter = 0.130, length = 0.075, r = {0, 0, 0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-124, 114}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder8(diameter = 0.130, length = 0.075, r = {0, 0, -0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-124, 160}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute4 annotation (
      Placement(visible = true, transformation(origin = {-169, 127}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.Fixed fixed(animation = false, r = {0.858, 0.065, 0.05}) annotation (
      Placement(visible = true, transformation(origin = {-234, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation(animation = true, r = {-0.153, 0.040, 0}) annotation (
      Placement(visible = true, transformation(origin = {-198, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation1(r = {-0.350, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-188, 22}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation2(r = {-0.25, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-190, -6}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute5 annotation (
      Placement(visible = true, transformation(origin = {-261, -125}, extent = {{-13, 13}, {13, -13}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation3(animation = false, r = {0.250, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-300, -126}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation4(r = {-0.192, -0.192, 0}) annotation (
      Placement(visible = true, transformation(origin = {-218, -124}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder9(diameter = 0.105, length = 0.075, r = {0, 0, -0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-132, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute6 annotation (
      Placement(visible = true, transformation(origin = {-181, -121}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder10(diameter = 0.105, length = 0.075, r = {0, 0, 0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-128, -148}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder3(diameter = 0.088, length = 0.075, r = {0, 0, -0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-122, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder11(diameter = 0.105, length = 0.075, r = {0, 0, 0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-134, -222}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute7 annotation (
      Placement(visible = true, transformation(origin = {-185, -203}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation5(animation = false, r = {0.281, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-310, -222}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder12(diameter = 0.105, length = 0.075, r = {0, 0, -0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-134, -182}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute8 annotation (
      Placement(visible = true, transformation(origin = {-276, -224}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation7(r = {-0.087, -0.029, 0}) annotation (
      Placement(visible = true, transformation(origin = {-288, 212}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder13(diameter = 0.105, length = 0.075, r = {0, 0, 0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-122, 194}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute9 annotation (
      Placement(visible = true, transformation(origin = {-175, 217}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder14(diameter = 0.105, length = 0.075, r = {0, 0, -0.0375}) annotation (
      Placement(visible = true, transformation(origin = {-120, 230}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation9(r = {-0.071, -0.120, 0}) annotation (
      Placement(visible = true, transformation(origin = {-208, 216}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute10 annotation (
      Placement(visible = true, transformation(origin = {-246, 218}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation6(r = {-0.192, -0.192, 0}) annotation (
      Placement(visible = true, transformation(origin = {-226, -202}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation8(animation = false, r = {-0.134, -0.134, 0}) annotation (
      Placement(visible = true, transformation(origin = {-288, -84}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation10(animation = false, r = {0.064, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-208, -48}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel springDamperParallel(c = 10000, d = 10000, fixedRotationAtFrame_a = false, fixedRotationAtFrame_b = false, s_unstretched = 0.15) annotation (
      Placement(visible = true, transformation(origin = {-264, -58}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation11(animation = false, r = {-0.134, -0.134, 0}) annotation (
      Placement(visible = true, transformation(origin = {-388, -182}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation12(animation = false, r = {0.098, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-376, -10}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel springDamperParallel1(c = 10000, d = 10000, fixedRotationAtFrame_a = false, fixedRotationAtFrame_b = false, s_unstretched = 0.15) annotation (
      Placement(visible = true, transformation(origin = {-401, -83}, extent = {{-11, -11}, {11, 11}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel springDamperParallel2(c = 10000, d = 10000, fixedRotationAtFrame_a = false, fixedRotationAtFrame_b = false, s_unstretched = 0.15) annotation (
      Placement(visible = true, transformation(origin = {-410, 56}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation13(animation = false, r = {-0.044, -0.065, 0}) annotation (
      Placement(visible = true, transformation(origin = {-250, 142}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation14(animation = false, r = {0.278, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-350, 30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    AGILEX_BUNKER_PRO.contact_sprocket contact_sprocket annotation (Placement(
          visible=true, transformation(
          origin={28,238},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    AGILEX_BUNKER_PRO.cinghiafinal cinghiafinal annotation (Placement(visible=
            true, transformation(
          origin={271,67},
          extent={{-67,-67},{67,67}},
          rotation=0)));
    Modelica.Mechanics.Rotational.Sources.Speed speed annotation (
      Placement(visible = true, transformation(origin = {-460, 268}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp(duration = 3, height = 5, startTime = 1) annotation (
      Placement(visible = true, transformation(origin = {-514, 268}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
    AGILEX_BUNKER_PRO.contact contact(R=0.0525) annotation (Placement(visible=
            true, transformation(
          origin={12,174},
          extent={{-40,-40},{40,40}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact1(R=0.065) annotation (Placement(visible=
            true, transformation(
          origin={6,78},
          extent={{-38,-38},{38,38}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact2(R=0.044) annotation (Placement(visible=
            true, transformation(
          origin={8,-10},
          extent={{-36,-36},{36,36}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact3(R=0.044) annotation (Placement(visible=
            true, transformation(
          origin={-10,-104},
          extent={{-44,-44},{44,44}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact4(R=0.044) annotation (Placement(visible=
            true, transformation(
          origin={-21,-193},
          extent={{-39,-39},{39,39}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact5(R=0.0525) annotation (Placement(visible=
            true, transformation(
          origin={-18,-288},
          extent={{-44,-44},{44,44}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact6(R=0.0525) annotation (Placement(visible=
            true, transformation(
          origin={-129,-347},
          extent={{-49,-49},{49,49}},
          rotation=0)));
  equation
    connect(revolute.frame_b, bodyCylinder.frame_a) annotation (
      Line(points = {{-176, 270}, {-140, 270}, {-140, 276}, {-130, 276}}));
    connect(bodyCylinder4.frame_a, revolute1.frame_b) annotation (
      Line(points = {{-134, 60}, {-134, 52.75}, {-138, 52.75}, {-138, 55.5}, {-156, 55.5}, {-156, 55}}, color = {95, 95, 95}));
    connect(bodyCylinder1.frame_a, revolute2.frame_b) annotation (
      Line(points = {{-134, 28}, {-146, 28}, {-146, 11}}, color = {95, 95, 95}));
    connect(bodyCylinder2.frame_a, revolute2.frame_b) annotation (
      Line(points = {{-136, -4}, {-136, -7.5}, {-146, -7.5}, {-146, 11}}, color = {95, 95, 95}));
    connect(bodyCylinder6.frame_a, revolute3.frame_b) annotation (
      Line(points = {{-142, -70}, {-142, -75.5}, {-162, -75.5}, {-162, -67}}, color = {95, 95, 95}));
    connect(bodyCylinder5.frame_a, revolute3.frame_b) annotation (
      Line(points = {{-140, -35}, {-140, -29.5}, {-166, -29.5}, {-166, -55.25}, {-162, -55.25}, {-162, -67}}, color = {95, 95, 95}));
    connect(bodyCylinder8.frame_a, revolute4.frame_b) annotation (
      Line(points = {{-134, 160}, {-134, 157.5}, {-158, 157.5}, {-158, 127}}, color = {95, 95, 95}));
    connect(bodyCylinder7.frame_a, revolute4.frame_b) annotation (
      Line(points = {{-134, 114}, {-134, 110.375}, {-158, 110.375}, {-158, 127}}, color = {95, 95, 95}));
    connect(fixed.frame_b, revolute4.frame_a) annotation (
      Line(points = {{-224, 46}, {-224, 99.5}, {-184, 99.5}, {-184, 100.25}, {-180, 100.25}, {-180, 127}}, color = {95, 95, 95}));
    connect(fixed.frame_b, fixedTranslation.frame_a) annotation (
      Line(points = {{-224, 46}, {-208, 46}}));
    connect(fixedTranslation.frame_b, revolute1.frame_a) annotation (
      Line(points = {{-188, 46}, {-188, 55}, {-178, 55}}));
    connect(fixedTranslation.frame_b, fixedTranslation1.frame_a) annotation (
      Line(points = {{-188, 46}, {-188, 32}}));
    connect(fixedTranslation1.frame_b, revolute2.frame_a) annotation (
      Line(points = {{-188, 12}, {-188, 11}, {-168, 11}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_b, revolute3.frame_a) annotation (
      Line(points = {{-190, -16}, {-191.25, -16}, {-191.25, -32}, {-192.5, -32}, {-192.5, -40}, {-184, -40}, {-184, -67}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_a, fixedTranslation1.frame_b) annotation (
      Line(points = {{-190, 4}, {-189, 4}, {-189, 12}, {-188, 12}}));
    connect(fixedTranslation3.frame_b, revolute5.frame_a) annotation (
      Line(points = {{-290, -126}, {-303, -126}, {-303, -125}, {-274, -125}}, color = {95, 95, 95}));
    connect(revolute5.frame_b, fixedTranslation4.frame_a) annotation (
      Line(points = {{-248, -125}, {-248, -137.5}, {-228, -137.5}, {-228, -124}}));
    connect(bodyCylinder10.frame_a, revolute6.frame_b) annotation (
      Line(points = {{-138, -148}, {-138, -151.813}, {-146, -151.813}, {-146, -149.625}, {-170, -149.625}, {-170, -121}}, color = {95, 95, 95}));
    connect(bodyCylinder9.frame_a, revolute6.frame_b) annotation (
      Line(points = {{-142, -110}, {-142, -106.75}, {-160, -106.75}, {-160, -105.5}, {-170, -105.5}, {-170, -121}}, color = {95, 95, 95}));
    connect(bodyCylinder3.frame_a, revolute1.frame_b) annotation (
      Line(points = {{-132, 86}, {-156, 86}, {-156, 55}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_b, fixedTranslation3.frame_a) annotation (
      Line(points = {{-190, -16}, {-190, -20}, {-310, -20}, {-310, -126}}, color = {95, 95, 95}));
    connect(bodyCylinder11.frame_a, revolute7.frame_b) annotation (
      Line(points = {{-144, -222}, {-144, -238.813}, {-174, -238.813}, {-174, -203}}, color = {95, 95, 95}));
    connect(bodyCylinder12.frame_a, revolute7.frame_b) annotation (
      Line(points = {{-146, -182}, {-146, -203.5}, {-174, -203.5}, {-174, -203}}, color = {95, 95, 95}));
    connect(fixedTranslation5.frame_b, revolute8.frame_a) annotation (
      Line(points = {{-300, -222}, {-300, -224}, {-286, -224}}, color = {95, 95, 95}));
    connect(fixedTranslation1.frame_b, fixedTranslation5.frame_a) annotation (
      Line(points = {{-188, 12}, {-320, 12}, {-320, -222}}, color = {95, 95, 95}));
    connect(fixed.frame_b, fixedTranslation7.frame_a) annotation (
      Line(points = {{-224, 46}, {-341.75, 46}, {-341.75, 176}, {-338.875, 176}, {-338.875, 212}, {-298, 212}}));
    connect(bodyCylinder13.frame_a, revolute9.frame_b) annotation (
      Line(points = {{-132, 194}, {-132, 185.188}, {-144, 185.188}, {-144, 184.375}, {-164, 184.375}, {-164, 217}}, color = {95, 95, 95}));
    connect(bodyCylinder14.frame_a, revolute9.frame_b) annotation (
      Line(points = {{-130, 230}, {-130, 240.5}, {-164, 240.5}, {-164, 217}}, color = {95, 95, 95}));
    connect(revolute9.frame_a, fixedTranslation9.frame_b) annotation (
      Line(points={{-186,217},{-198,217},{-198,216}},        color = {95, 95, 95}));
    connect(revolute10.frame_b, fixedTranslation9.frame_a) annotation (
      Line(points = {{-236, 218}, {-223, 218}, {-223, 216}, {-218, 216}}));
    connect(fixedTranslation7.frame_b, revolute10.frame_a) annotation (
      Line(points = {{-278, 212}, {-256, 212}, {-256, 218}}, color = {95, 95, 95}));
    connect(revolute8.frame_b, fixedTranslation6.frame_a) annotation (
      Line(points = {{-266, -224}, {-240.5, -224}, {-240.5, -202}, {-236, -202}}));
    connect(revolute7.frame_a, fixedTranslation6.frame_b) annotation (
      Line(points = {{-196, -203}, {-216, -203}, {-216, -202}}, color = {95, 95, 95}));
    connect(revolute6.frame_a, fixedTranslation4.frame_b) annotation (
      Line(points = {{-192, -121}, {-208, -121}, {-208, -124}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_b, fixedTranslation10.frame_a) annotation (
      Line(points = {{-190, -16}, {-190, -32}, {-198, -32}, {-198, -48}}, color = {95, 95, 95}));
    connect(fixedTranslation8.frame_a, revolute5.frame_b) annotation (
      Line(points={{-298,-84},{-298,-114},{-248,-114},{-248,-125}}));
    connect(springDamperParallel.frame_b, fixedTranslation8.frame_b) annotation (
      Line(points = {{-264, -68}, {-252, -68}, {-252, -84}, {-278, -84}}, color = {95, 95, 95}));
    connect(springDamperParallel.frame_a, fixedTranslation10.frame_b) annotation (
      Line(points = {{-264, -48}, {-218, -48}}, color = {95, 95, 95}));
    connect(fixedTranslation12.frame_a, fixedTranslation1.frame_b) annotation (
      Line(points = {{-366, -10}, {-220, -10}, {-220, 12}, {-188, 12}}));
    connect(fixedTranslation12.frame_b, springDamperParallel1.frame_a) annotation (
      Line(points = {{-386, -10}, {-401, -10}, {-401, -72}}, color = {95, 95, 95}));
    connect(fixedTranslation11.frame_a, revolute8.frame_b) annotation (
      Line(points = {{-378, -182}, {-266, -182}, {-266, -224}}, color = {95, 95, 95}));
    connect(fixedTranslation11.frame_b, springDamperParallel1.frame_b) annotation (
      Line(points = {{-398, -182}, {-401, -182}, {-401, -94}}, color = {95, 95, 95}));
    connect(fixedTranslation14.frame_a, fixedTranslation1.frame_b) annotation (
      Line(points = {{-340, 30}, {-188, 30}, {-188, 12}}, color = {95, 95, 95}));
    connect(fixedTranslation14.frame_b, springDamperParallel2.frame_a) annotation (
      Line(points = {{-360, 30}, {-410, 30}, {-410, 46}}, color = {95, 95, 95}));
    connect(revolute10.frame_b, fixedTranslation13.frame_a) annotation (
      Line(points = {{-236, 218}, {-224, 218}, {-224, 152}, {-250, 152}}));
    connect(fixedTranslation13.frame_b, springDamperParallel2.frame_b) annotation (
      Line(points = {{-250, 132}, {-410, 132}, {-410, 66}}));
    connect(bodyCylinder.frame_b, contact_sprocket.frame_b) annotation (
      Line(points={{-110,276},{16.2,276},{16.2,238.4}}));
    connect(speed.flange, revolute.axis) annotation (
      Line(points = {{-450, 268}, {-254, 268}, {-254, 290}, {-186, 290}, {-186, 280}}));
    connect(ramp.y, speed.w_ref) annotation (
      Line(points={{-505.2,268},{-472,268}},    color = {0, 0, 127}));
    connect(revolute9.frame_b, contact.frame_b) annotation (
      Line(points={{-164,217},{-72,217},{-72,175.6},{-35.2,175.6}}));
    connect(contact.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{52,174.8},{204,174.8},{204,65.66}}, color = {95, 95, 95}, thickness = 0.5));
    connect(revolute4.frame_b, contact1.frame_b) annotation (
      Line(points={{-158,127},{-60,127},{-60,79.52},{-38.84,79.52}}, color = {95, 95, 95}));
    connect(contact1.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{44,78.76},{204,78.76},{204,65.66}},
                                                      thickness = 0.5));
    connect(revolute1.frame_b, contact2.frame_b) annotation (
      Line(points={{-156,55},{-58,55},{-58,-8.56},{-34.48,-8.56}}));
    connect(contact2.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{44,-9.28},{204,-9.28},{204,65.66}}, thickness = 0.5));
    connect(revolute2.frame_b, contact3.frame_b) annotation (
      Line(points={{-146,11},{-76,11},{-76,-102.24},{-61.92,-102.24}}, color = {95, 95, 95}));
    connect(contact3.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{34,-103.12},{204,-103.12},{204,65.66}},
                                                          color = {95, 95, 95}, thickness = 0.5));
    connect(revolute3.frame_b, contact4.frame_b) annotation (
      Line(points={{-162,-67},{-82,-67},{-82,-191.44},{-67.02,-191.44}}));
    connect(contact4.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{18,-192.22},{204,-192.22},{204,65.66}},
                                                          thickness = 0.5));
    connect(revolute6.frame_b, contact5.frame_b) annotation (
      Line(points={{-170,-121},{-106,-121},{-106,-286.24},{-69.92,-286.24}}, color = {95, 95, 95}));
    connect(contact5.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{26,-287.12},{204,-287.12},{204,65.66}},
                                                          color = {95, 95, 95}, thickness = 0.5));
    connect(revolute7.frame_b, contact6.frame_b) annotation (
      Line(points={{-174,-203},{-224,-203},{-224,-345.04},{-186.82,-345.04}}));
    connect(contact6.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{-80,-346.02},{204,-346.02},{204,65.66}},
                                                           thickness = 0.5));
    connect(contact_sprocket.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{38,238.2},{124,238.2},{124,65.66},{204,65.66}},thickness = 0.5));
    connect(revolute.frame_a, world.frame_b) annotation (
      Line(points = {{-196, 270}, {-196, 264.25}, {-228, 264.25}, {-228, 236.5}, {-266, 236.5}, {-266, 241}, {-358, 241}, {-358, 186}}));
    annotation (
      Icon(coordinateSystem(grid = {2, 0})),
      experiment(StartTime = 0, StopTime = 5, Tolerance = 0.001, Interval = 0.01));
  end trackfinal;

  model cinghiafinal
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b[25] annotation (
      Placement(visible = true, transformation(origin = {-100, -3}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-100, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a[25] annotation (
      Placement(visible = true, transformation(origin = {100, -9}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {100, -8}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    AGILEX_BUNKER_PRO.element link[25](
      startX={0.583867,0.671173,0.755162,0.836694,0.916869,0.930544,0.849205,
          0.757759,0.666175,0.574575,0.482967,0.39136,0.299764,0.20819,0.116621,
          0.030362,-0.0483488,-0.0935338,-0.0486413,0.0385476,0.126538,0.217664,
          0.309245,0.400848,0.492266},
      startY={-0.138365,-0.11062,-0.0740411,-0.0322689,0.0120491,0.102622,
          0.144752,0.150211,0.152326,0.151047,0.150645,0.151221,0.152791,
          0.155327,0.158032,0.127189,0.0803199,0.000634477,-0.079213,-0.107326,
          -0.132819,-0.14222,-0.144497,-0.145531,-0.139609},
      startZ=fill(0.05, 25),
      thZ={0.30771657,0.410737918,0.473475722,0.504965371,1.420869243,
          2.663644401,3.081961127,3.118500936,3.155556,3.145992237,3.135304458,
          3.124452162,3.113914473,3.112068688,3.484972025,3.678667741,
          4.196513444,5.224501201,5.971263171,6.001175396,6.18037405,
          6.258324479,6.271896838,6.347869744,6.296755677}) annotation (
        Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));
  equation
    for n in 1:24 loop
      connect(link[n].frame_b, link[n + 1].frame_a);
    end for;
    connect(link[25].frame_b, link[1].frame_a);
    for n in 1:25 loop
      connect(link[n].frame_b1, frame_b[n]);
      connect(link[n].frame_a1, frame_a[n]);
    end for;
    annotation (
      Icon(coordinateSystem(grid = {2, 0})),
      Diagram(coordinateSystem(grid = {2, 3})));
  end cinghiafinal;

  model contact_force_idler_new
   Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a annotation (
      Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b annotation (
      Placement(visible = true, transformation(origin = {100, 0}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {100, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Units.SI.Length DeltaD[3];
    Modelica.Units.SI.Velocity dDeltaD[3];
    Modelica.Units.SI.Distance dist;
    Modelica.Mechanics.MultiBody.Frames.Orientation Rrel;
    Modelica.Units.SI.Length DeltaDrel[3];
    Modelica.Units.SI.Velocity dDeltaDrel[3];
    Modelica.Units.SI.AngularVelocity rw[3];
    Modelica.Units.SI.Angle rth[2], rth1[3];
    Modelica.Units.SI.Velocity vrel;
    Modelica.Units.SI.Force Fn, Ft, Fl;
    Modelica.Units.SI.Torque T2, T3, T[3];
    Real smooth1;
    Real E[3,3];
    parameter Modelica.Units.SI.Length R = 0.075;
    parameter Modelica.Units.SI.Length len = (0.0586*38)/80;
    parameter Modelica.Units.SI.Length lar = 0.15;
    parameter Modelica.Units.SI.Length L = 0.048;
    parameter Modelica.Units.SI.ModulusOfElasticity Eal = 70e9;
    parameter Modelica.Units.SI.ModulusOfElasticity Er = 1e7;
    parameter Modelica.Units.SI.PoissonNumber val = 0.33;
    parameter Modelica.Units.SI.PoissonNumber vr = 0.5;
    final parameter Modelica.Units.SI.ModulusOfElasticity Estar = 1/((1-vr^2)/Er + (1-val^2)/Eal);
    final parameter Real KR = (Modelica.Constants.pi/4)*Estar*L;//7*10 ^ 5;
    final parameter Real DR = KR/1000;//10 ^ 2;
    //parameter Real KT = 7*10 ^ 5;
    //parameter Real DT = 10 ^ 2;
    final parameter Real KL = KR;//7*10 ^ 5;
    final parameter Real DL = KL/1000;//10 ^ 2;
    parameter Real Krx = KR*(lar/2 - L/4)^2;//7*10 ^ 5;
    parameter Real Drx = DR*(lar/2 - L/4)^2;//10 ^ 2;
    parameter Real Kry = Krx;//7*10 ^ 5;
    parameter Real Dry = Kry/1000;//10 ^ 2;
    //parameter Modelica.Units.SI.Length lar = 0.105;
    //parameter Modelica.Units.SI.Length s = 0.0125;
    //Friction coefficient function

    function mu
      input Modelica.Units.SI.Velocity v_rel;
      input Modelica.Units.SI.CoefficientOfFriction mue_s = 0.4;
      input Modelica.Units.SI.CoefficientOfFriction mue_k = 0.3;
      input Modelica.Units.SI.Velocity v_e1 = 0.01;
      input Modelica.Units.SI.Velocity v_e2 = 0.1;
      input Real k_v(unit = "s2/m") = 0;
      output Real m;
    protected
      Real gamma1;
      Real gamma2;
      Real gamma3;
      Real gamma4;
      Real gamma5;
      Real gamma6;
    algorithm
      gamma4 := mue_k - k_v * v_e2;
      gamma1 := (mue_s - tanh(2 * v_e1 / v_e2) * gamma4 - k_v * 2 * v_e1) / (tanh(2 * v_e2 / v_e1) - tanh(2 * v_e1 / v_e2)) - gamma4;
      gamma2 := 2 / v_e1;
      gamma3 := 3 / v_e2;
      gamma5 := gamma2;
      gamma6 := k_v;
      m := gamma1 * (tanh(gamma2 * abs(v_rel)) - tanh(gamma3 * abs(v_rel))) + gamma4 * tanh(gamma5 * abs(v_rel)) + gamma6 * abs(v_rel);
    end mu;
  equation
    //Relative rotation between frame_a and frame_b
    Rrel = Modelica.Mechanics.MultiBody.Frames.relativeRotation(frame_a.R, frame_b.R);
    //Distance between the center of frame_a and the center of frame_b
    dist = Modelica.Math.Vectors.length(DeltaD);
    //Position vector of the center of frame_a from the center of frame_b resolved in world reference
    DeltaD = frame_a.r_0 - frame_b.r_0;
    //Relative velocity of frame_a wrt frame_b resolved in world reference
    dDeltaD = der(frame_a.r_0) - der(frame_b.r_0);
    //Position vector of the center of frame_b from the center of frame_a resolved in frame_a reference
    DeltaDrel = Modelica.Mechanics.MultiBody.Frames.resolve2(frame_a.R, -DeltaD);
    //Relative velocity of frame_b wrt frame_a resolved in frame_a reference
    dDeltaDrel = Modelica.Mechanics.MultiBody.Frames.resolve2(frame_a.R, -dDeltaD);
    //Relative angular velocities and angles of frame_b wrt frame_a resolved in frame_a
    rw = -Modelica.Mechanics.MultiBody.Frames.resolve1(Rrel, {Rrel.w[1], Rrel.w[2], 0});
    E = [1, (Modelica.Math.sin(rth[1])*Modelica.Math.sin(rth[2]))/Modelica.Math.cos(rth[2]), -(Modelica.Math.cos(rth[1])*Modelica.Math.sin(rth[2]))/Modelica.Math.cos(rth[2]);
       0, Modelica.Math.cos(rth[1]), Modelica.Math.sin(rth[1]);
       0, -Modelica.Math.sin(rth[1])/Modelica.Math.cos(rth[2]), Modelica.Math.cos(rth[1])/Modelica.Math.cos(rth[2])];
    der(rth) = {E[1,1:3]*rw, E[2,1:3]*rw};
    rth1 = - Modelica.Mechanics.MultiBody.Frames.axesRotationsAngles(Rrel, {1,2,3});

    when Modelica.Math.Vectors.length(rth - rth1[1:2]) > 1e-5 then
    reinit(rth, rth1[1:2]);
    end when;
    //Slip velocity of the cylinder rolling on the surface of the element
    vrel = dDeltaDrel[1] + frame_b.R.w[3] * R;
    //Components of the total force acting on the element
    Fn = smooth1 * smooth(0, max(0, -(KR * (DeltaDrel[2] - R) + DR * dDeltaDrel[2])));
    Ft = -mu(vrel) * Fn * sign(vrel);
    Fl = smooth1 * (-(DL * dDeltaDrel[3] + KL * DeltaDrel[3]));
    frame_a.f = {Ft, Fn, Fl};
    //Components of the total torque acting on the element
    T = transpose(E)*{smooth1 * (Drx * der(rth[1]) + Krx * rth[1]), smooth1 * (Dry * der(rth[2]) + Kry * rth[2]), 0};
    T2 = T[2] - DeltaDrel[1] * Fl;
    T3 = DeltaDrel[1] * Fn + T[3];
    frame_a.t = {T[1], T2, T3};
    //Equilibrium of forces and torques acting on the cylinder
    frame_b.f + Modelica.Mechanics.MultiBody.Frames.resolve2(Rrel, frame_a.f) = {0, 0, 0};
    frame_b.t + Modelica.Mechanics.MultiBody.Frames.resolve2(Rrel, {T[1] - R * Fl, T[2], R * Ft}) = {0, 0, 0};
  algorithm
    //Computation of the smooth function
    //smooth := ((Modelica.Math.atan(-100000 * (DeltaDrel[2] - R)) + Modelica.Constants.pi / 2) / Modelica.Constants.pi) * ((Modelica.Math.atan(-100000 * (abs(DeltaDrel[1]) - len / 2)) + Modelica.Constants.pi / 2) / Modelica.Constants.pi);
    smooth1 := ((Modelica.Math.tanh(-1000 * (DeltaDrel[2] - R)) + 1) / 2)* ((Modelica.Math.tanh(-1000 * (abs(DeltaDrel[1]) - len / 2)) + 1) / 2);
  initial equation
    rth =  rth1[1:2];
    annotation (
      Icon(coordinateSystem(grid = {2, 0})),
      experiment(StartTime = 0, StopTime = 0.3, Tolerance = 1e-6, Interval = 0.0006));
  end contact_force_idler_new;

  model contact_force_new
    Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a annotation (
      Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-100, 0}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b annotation (
      Placement(visible = true, transformation(origin = {100, 0}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {100, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Units.SI.Length DeltaD[3];
    Modelica.Units.SI.Velocity dDeltaD[3];
    Modelica.Units.SI.Distance dist;
    Modelica.Mechanics.MultiBody.Frames.Orientation Rrel;
    Modelica.Units.SI.Length DeltaDrel[3];
    Modelica.Units.SI.Velocity dDeltaDrel[3];
    Modelica.Units.SI.AngularVelocity rw[3];
    Modelica.Units.SI.Angle rth[2], rth1[3];
    Real smooth1;
    Real E[3,3];
    //, smooth1;
    Modelica.Units.SI.Force Fn, Ft, Fl;
    Modelica.Units.SI.Torque T2, T3,T[3];
    Modelica.Units.SI.Velocity vrel;
    //Modelica.Units.SI.Angle relangles[3];
    //Modelica.Mechanics.MultiBody.Frames.Orientation Rn;
    //Modelica.Mechanics.MultiBody.Frames.Orientation Rrel1;
    parameter Modelica.Units.SI.Length R = 0.08;
    parameter Modelica.Units.SI.Length len = (0.0586*38)/80;
    parameter Modelica.Units.SI.Length lar = 0.105;
    parameter Modelica.Units.SI.Length L = lar;
    parameter Modelica.Units.SI.ModulusOfElasticity Eal = 70e9;
    parameter Modelica.Units.SI.ModulusOfElasticity Er = 1e7;
    parameter Modelica.Units.SI.PoissonNumber val = 0.33;
    parameter Modelica.Units.SI.PoissonNumber vr = 0.5;
    final parameter Modelica.Units.SI.ModulusOfElasticity Estar = 1/((1-vr^2)/Er + (1-val^2)/Eal);
    final parameter Real KR = (Modelica.Constants.pi/4)*Estar*L;//7*10 ^ 5;
    final parameter Real DR = KR/1000;//10 ^ 2;
    //parameter Real KT = 7*10 ^ 5;
    //parameter Real DT = 10 ^ 2;
    final parameter Real KL = KR;//7*10 ^ 5;
    final parameter Real DL = KL/1000;//10 ^ 2;
    final parameter Real Krx = KR*(lar/2 - L/4)^2;//7*10 ^ 5;
    final parameter Real Drx = DR*(lar/2 - L/4)^2;//10 ^ 2;
    final parameter Real Kry = Krx;//7*10 ^ 5;
    final parameter Real Dry = Kry/1000;//10 ^ 2;
    //parameter Real KR = 7*10 ^ 5;
    //parameter Real DR = 10 ^ 2;
    //parameter Real KT = 10 ^ 6;
    parameter Real DT = 7*10 ^ 5;
    //parameter Real KL = 7*10 ^ 5;
    //parameter Real DL = 10 ^ 2;
    //parameter Real Krx = 7*10 ^ 5;
    //parameter Real Drx = 10 ^ 2;
    //parameter Real Kry = 7*10 ^ 5;
    //parameter Real Dry = 10 ^ 2;

    //parameter Modelica.Units.SI.Length lar = 0.105;
  equation
    //Relative rotation between frame_a and frame_b
    Rrel = Modelica.Mechanics.MultiBody.Frames.relativeRotation(frame_a.R, frame_b.R);
    //Relative rotation of normal reference wrt world reference
    //Rn = Modelica.Mechanics.MultiBody.Frames.from_nxy(cross(-Modelica.Math.Vectors.normalize(DeltaD), {0, 0, 1}), -Modelica.Math.Vectors.normalize(DeltaD));
    //Relative rotation of normal reference wrt frame_a
    //Rrel1 = Modelica.Mechanics.MultiBody.Frames.relativeRotation(frame_a.R, Rn);
    //Distance between the center of frame_a and the center of frame_b
    dist = Modelica.Math.Vectors.length(DeltaD);
    //Position vector of the center of frame_a from the center of frame_b resolved in world reference
    DeltaD = frame_a.r_0 - frame_b.r_0;
    //Relative velocity of frame_a wrt frame_b resolved in world reference
    dDeltaD = der(frame_a.r_0) - der(frame_b.r_0);
    //Velocità relativa del frame_a rispetto al frame_b riportato frame_n
    DeltaDrel = Modelica.Mechanics.MultiBody.Frames.resolve2(frame_a.R, -DeltaD);
    //Relative velocity of frame_b wrt frame_a resolved in frame_a reference
    dDeltaDrel = Modelica.Mechanics.MultiBody.Frames.resolve2(frame_a.R, -dDeltaD);
    //Relative angular velocities and angles of frame_b wrt frame_a resolved in frame_a
    rw = -Modelica.Mechanics.MultiBody.Frames.resolve1(Rrel, {Rrel.w[1], Rrel.w[2], 0});
    E = [1, (Modelica.Math.sin(rth[1])*Modelica.Math.sin(rth[2]))/Modelica.Math.cos(rth[2]), -(Modelica.Math.cos(rth[1])*Modelica.Math.sin(rth[2]))/Modelica.Math.cos(rth[2]);
       0, Modelica.Math.cos(rth[1]), Modelica.Math.sin(rth[1]);
       0, -Modelica.Math.sin(rth[1])/Modelica.Math.cos(rth[2]), Modelica.Math.cos(rth[1])/Modelica.Math.cos(rth[2])];

    der(rth) = {E[1,1:3]*rw, E[2,1:3]*rw};//E*rw;
    rth1 = - Modelica.Mechanics.MultiBody.Frames.axesRotationsAngles(Rrel, {1,2,3});

    when Modelica.Math.Vectors.length(rth - rth1[1:2]) > 1e-5 then
    reinit(rth, rth1[1:2]);
    end when;
    //Relative angles of frame_b wrt frame_a resolved in frame_a
    //relangles = -Modelica.Mechanics.MultiBody.Frames.axesRotationsAngles(Rrel1, {1, 2, 3});
    //Components of the total force acting on the element
    //if noEvent(dist - R < 0) then -(DT * (dDeltaDrel[1] + frame_b.R.w[3] * R) + KT * DeltaDrel[1]) else 0
    //smooth1 * (-(DT * (dDeltaDrel[1] + frame_b.R.w[3] * R) + KT * DeltaDrel[1]))
    //Ft = smooth1 * (-(DT * (dDeltaDrel[1] + frame_b.R.w[3] * R) + KT * DeltaDrel[1]));
    //Slip velocity of the cylinder rolling on the surface of the element
    vrel = dDeltaDrel[1] + frame_b.R.w[3] * R;
    Ft = smooth1 * (-DT * vrel);
    Fn = smooth1 * smooth(0, max(0, -(KR * (DeltaDrel[2] - R) + DR * dDeltaDrel[2])));
    Fl = smooth1 * (-(DL * dDeltaDrel[3] + KL * DeltaDrel[3]));
    frame_a.f = {Ft, Fn, Fl};
    //Components of the total torque acting on the element
    T = transpose(E)*{smooth1 * (Drx * der(rth[1]) + Krx * rth[1]), smooth1 * (Dry * der(rth[2]) + Kry * rth[2]), 0};
    T2 = T[2] - DeltaDrel[1] * Fl;
    T3 = DeltaDrel[1] * Fn + T[3];
    frame_a.t = {T[1], T2, T3};
    //Equilibrium of forces and torques acting on the cylinder
    frame_b.f + Modelica.Mechanics.MultiBody.Frames.resolve2(Rrel, frame_a.f) = {0, 0, 0};
    frame_b.t + Modelica.Mechanics.MultiBody.Frames.resolve2(Rrel, {T[1] - R * Fl, T[2], R * Ft}) = {0, 0, 0};
    //frame_b.t + Modelica.Mechanics.MultiBody.Frames.resolve2(Rrel, {frame_a.t[1] + R * frame_a.f[3], smooth * (Dr * rw[2] + Kr * rth[2]), R * frame_a.f[1]}) = {0, 0, 0};
  algorithm
    //Computation of the smooth function
    //smooth := ((Modelica.Math.atan(-100000 * (DeltaDrel[2] - R)) + Modelica.Constants.pi / 2) / Modelica.Constants.pi) * ((Modelica.Math.atan(-100000 * (abs(DeltaDrel[1]) - len / 2)) + Modelica.Constants.pi / 2) / Modelica.Constants.pi);
    smooth1 := ((Modelica.Math.tanh(-1000 * (DeltaDrel[2] - R)) + 1) / 2)* ((Modelica.Math.tanh(-1000 * (abs(DeltaDrel[1]) - len / 2)) + 1) / 2);
    //smooth1 := ((Modelica.Math.tanh(-500 * ((dist - R)-0.01)) + 1) / 2);
    //smooth1 := (Modelica.Math.tanh(-10000 * (dist - R)) + 1) / 2;
    //smooth1 := (Modelica.Math.atan(-30000 * ((dist - R)-0.01)) + Modelica.Constants.pi / 2) / Modelica.Constants.pi;
  initial equation
    rth =  rth1[1:2];
    annotation (
      Icon(coordinateSystem(grid = {2, 0})),
      experiment(StartTime = 0, StopTime = 0.3, Tolerance = 1e-6, Interval = 0.0006));
  end contact_force_new;

  model trackterrain
    extends Modelica.Icons.Example;
    inner Modelica.Mechanics.MultiBody.World world annotation (
      Placement(visible = true, transformation(origin = {-506, 174}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder(diameter = 0.16, length = 0.1, lengthDirection = {0, 0, 1}, r = {0, 0, 0.05}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-120, 276}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute(animation = false, useAxisFlange = true) annotation (
      Placement(visible = true, transformation(origin = {-186, 270}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder4(diameter = 0.088, length = 0.075, r = {0, 0, 0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-124, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute1 annotation (
      Placement(visible = true, transformation(origin = {-167, 55}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder1(diameter = 0.088, length = 0.075, r = {0, 0, -0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-124, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder2(diameter = 0.088, length = 0.075, r = {0, 0, 0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-126, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute2 annotation (
      Placement(visible = true, transformation(origin = {-157, 11}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder5(diameter = 0.088, length = 0.075, r = {0, 0, -0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-129, -35}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder6(diameter = 0.088, length = 0.075, r = {0, 0, 0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-132, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute3 annotation (
      Placement(visible = true, transformation(origin = {-173, -67}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder7(diameter = 0.130, length = 0.075, r = {0, 0, 0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-124, 114}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder8(diameter = 0.130, length = 0.075, r = {0, 0, -0.0375}, r_0(each fixed = true, start = {0.858, 0.065 + 0.2, 0.05}), useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-124, 160}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute4 annotation (
      Placement(visible = true, transformation(origin = {-169, 127}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation(animation = true, r = {-0.153, 0.040, 0}) annotation (
      Placement(visible = true, transformation(origin = {-198, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation1(r = {-0.350, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-188, 22}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation2(r = {-0.25, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-190, -6}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute5 annotation (
      Placement(visible = true, transformation(origin = {-261, -125}, extent = {{-13, 13}, {13, -13}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation3(animation = false, r = {0.250, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-300, -126}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation4(r = {-0.192, -0.192, 0}) annotation (
      Placement(visible = true, transformation(origin = {-218, -124}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder9(diameter = 0.105, length = 0.075, r = {0, 0, -0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-132, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute6 annotation (
      Placement(visible = true, transformation(origin = {-181, -121}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder10(diameter = 0.105, length = 0.075, r = {0, 0, 0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-128, -148}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder3(diameter = 0.088, length = 0.075, r = {0, 0, -0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-122, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder11(diameter = 0.105, length = 0.075, r = {0, 0, 0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-134, -222}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute7 annotation (
      Placement(visible = true, transformation(origin = {-185, -203}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation5(animation = false, r = {0.281, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-310, -222}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder12(diameter = 0.105, length = 0.075, r = {0, 0, -0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-134, -182}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute8 annotation (
      Placement(visible = true, transformation(origin = {-276, -224}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation7(r = {-0.087, -0.029, 0}) annotation (
      Placement(visible = true, transformation(origin = {-288, 212}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder13(diameter = 0.105, length = 0.075, r = {0, 0, 0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-122, 194}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute9 annotation (
      Placement(visible = true, transformation(origin = {-175, 217}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder14(diameter = 0.105, length = 0.075, r = {0, 0, -0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-120, 230}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation9(r = {-0.071, -0.120, 0}) annotation (
      Placement(visible = true, transformation(origin = {-208, 216}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute10 annotation (
      Placement(visible = true, transformation(origin = {-246, 218}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation6(r = {-0.192, -0.192, 0}) annotation (
      Placement(visible = true, transformation(origin = {-226, -202}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation8(animation = false, r = {-0.134, -0.134, 0}) annotation (
      Placement(visible = true, transformation(origin = {-288, -84}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation10(animation = false, r = {0.064, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-208, -48}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel springDamperParallel(c = 100000, d = 100000, fixedRotationAtFrame_a = false, fixedRotationAtFrame_b = false, s_unstretched = 0.15) annotation (
      Placement(visible = true, transformation(origin = {-264, -58}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation11(animation = false, r = {-0.134, -0.134, 0}) annotation (
      Placement(visible = true, transformation(origin = {-388, -182}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation12(animation = false, r = {0.098, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-376, -10}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel springDamperParallel1(c = 100000, d = 100000, fixedRotationAtFrame_a = false, fixedRotationAtFrame_b = false, s_unstretched = 0.15) annotation (
      Placement(visible = true, transformation(origin = {-401, -83}, extent = {{-11, -11}, {11, 11}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel springDamperParallel2(c = 100000, d = 100000, fixedRotationAtFrame_a = false, fixedRotationAtFrame_b = false, s_unstretched = 0.15) annotation (
      Placement(visible = true, transformation(origin = {-410, 56}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation13(animation = false, r = {-0.044, -0.065, 0}) annotation (
      Placement(visible = true, transformation(origin = {-250, 142}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation14(animation = false, r = {0.278, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-350, 30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    AGILEX_BUNKER_PRO.contact_sprocket contact_sprocket annotation (Placement(
          visible=true, transformation(
          origin={28,238},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    AGILEX_BUNKER_PRO.cinghiaterrain cinghiafinal(start={0,0.20,0}) annotation
      (Placement(visible=true, transformation(
          origin={269,67},
          extent={{-67,-67},{67,67}},
          rotation=0)));
    Modelica.Mechanics.Rotational.Sources.Speed speed annotation (
      Placement(visible = true, transformation(origin = {-460, 268}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp(duration = 3, height = 3, startTime = 1) annotation (
      Placement(visible = true, transformation(origin = {-514, 268}, extent = {{-8, -8}, {8, 8}}, rotation = 0)));
    AGILEX_BUNKER_PRO.contact contact(R=0.0525) annotation (Placement(visible=
            true, transformation(
          origin={12,174},
          extent={{-40,-40},{40,40}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact1(R=0.065) annotation (Placement(visible=
            true, transformation(
          origin={6,78},
          extent={{-38,-38},{38,38}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact2(R=0.044) annotation (Placement(visible=
            true, transformation(
          origin={8,-10},
          extent={{-36,-36},{36,36}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact3(R=0.044) annotation (Placement(visible=
            true, transformation(
          origin={-10,-104},
          extent={{-44,-44},{44,44}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact4(R=0.044) annotation (Placement(visible=
            true, transformation(
          origin={-21,-193},
          extent={{-39,-39},{39,39}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact5(R=0.0525) annotation (Placement(visible=
            true, transformation(
          origin={-18,-288},
          extent={{-44,-44},{44,44}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact6(R=0.0525) annotation (Placement(visible=
            true, transformation(
          origin={-129,-347},
          extent={{-49,-49},{49,49}},
          rotation=0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation15(r = {-0.105, -0.105, -0.05}) annotation (
      Placement(visible = true, transformation(origin = {-562, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    AGILEX_BUNKER_PRO.contact_t contact_t annotation (Placement(visible=true,
          transformation(
          origin={562,60},
          extent={{-50,-50},{50,50}},
          rotation=0)));
    inner AGILEX_BUNKER_PRO.Surface Surface(LengthX=20, LengthZ=20) annotation
      (Placement(visible=true, transformation(
          origin={398,272},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    Modelica.Mechanics.MultiBody.Forces.Force force annotation (
      Placement(visible = true, transformation(origin = {-206, 326}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp1(duration = 2, height = 10000.0, startTime = 2) annotation (
      Placement(visible = true, transformation(origin = {-830, 166}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.RealExpression realExpression annotation (
      Placement(visible = true, transformation(origin = {-722, 302}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.RealExpression realExpression1 annotation (
      Placement(visible = true, transformation(origin = {-614, 332}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(revolute.frame_b, bodyCylinder.frame_a) annotation (
      Line(points = {{-176, 270}, {-140, 270}, {-140, 276}, {-130, 276}}));
    connect(bodyCylinder4.frame_a, revolute1.frame_b) annotation (
      Line(points = {{-134, 60}, {-134, 52.75}, {-138, 52.75}, {-138, 55.5}, {-156, 55.5}, {-156, 55}}, color = {95, 95, 95}));
    connect(bodyCylinder1.frame_a, revolute2.frame_b) annotation (
      Line(points = {{-134, 28}, {-146, 28}, {-146, 11}}, color = {95, 95, 95}));
    connect(bodyCylinder2.frame_a, revolute2.frame_b) annotation (
      Line(points = {{-136, -4}, {-136, -7.5}, {-146, -7.5}, {-146, 11}}, color = {95, 95, 95}));
    connect(bodyCylinder6.frame_a, revolute3.frame_b) annotation (
      Line(points = {{-142, -70}, {-142, -75.5}, {-162, -75.5}, {-162, -67}}, color = {95, 95, 95}));
    connect(bodyCylinder5.frame_a, revolute3.frame_b) annotation (
      Line(points = {{-140, -35}, {-140, -29.5}, {-166, -29.5}, {-166, -55.25}, {-162, -55.25}, {-162, -67}}, color = {95, 95, 95}));
    connect(bodyCylinder8.frame_a, revolute4.frame_b) annotation (
      Line(points = {{-134, 160}, {-134, 157.5}, {-158, 157.5}, {-158, 127}}, color = {95, 95, 95}));
    connect(bodyCylinder7.frame_a, revolute4.frame_b) annotation (
      Line(points = {{-134, 114}, {-134, 110.375}, {-158, 110.375}, {-158, 127}}, color = {95, 95, 95}));
    connect(fixedTranslation.frame_b, revolute1.frame_a) annotation (
      Line(points = {{-188, 46}, {-188, 55}, {-178, 55}}));
    connect(fixedTranslation.frame_b, fixedTranslation1.frame_a) annotation (
      Line(points = {{-188, 46}, {-188, 32}}));
    connect(fixedTranslation1.frame_b, revolute2.frame_a) annotation (
      Line(points = {{-188, 12}, {-188, 11}, {-168, 11}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_b, revolute3.frame_a) annotation (
      Line(points = {{-190, -16}, {-191.25, -16}, {-191.25, -32}, {-192.5, -32}, {-192.5, -40}, {-184, -40}, {-184, -67}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_a, fixedTranslation1.frame_b) annotation (
      Line(points = {{-190, 4}, {-189, 4}, {-189, 12}, {-188, 12}}));
    connect(fixedTranslation3.frame_b, revolute5.frame_a) annotation (
      Line(points = {{-290, -126}, {-303, -126}, {-303, -125}, {-274, -125}}, color = {95, 95, 95}));
    connect(revolute5.frame_b, fixedTranslation4.frame_a) annotation (
      Line(points = {{-248, -125}, {-248, -137.5}, {-228, -137.5}, {-228, -124}}));
    connect(bodyCylinder10.frame_a, revolute6.frame_b) annotation (
      Line(points = {{-138, -148}, {-138, -151.813}, {-146, -151.813}, {-146, -149.625}, {-170, -149.625}, {-170, -121}}, color = {95, 95, 95}));
    connect(bodyCylinder9.frame_a, revolute6.frame_b) annotation (
      Line(points = {{-142, -110}, {-142, -106.75}, {-160, -106.75}, {-160, -105.5}, {-170, -105.5}, {-170, -121}}, color = {95, 95, 95}));
    connect(bodyCylinder3.frame_a, revolute1.frame_b) annotation (
      Line(points = {{-132, 86}, {-156, 86}, {-156, 55}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_b, fixedTranslation3.frame_a) annotation (
      Line(points = {{-190, -16}, {-190, -20}, {-310, -20}, {-310, -126}}, color = {95, 95, 95}));
    connect(bodyCylinder11.frame_a, revolute7.frame_b) annotation (
      Line(points = {{-144, -222}, {-144, -238.813}, {-174, -238.813}, {-174, -203}}, color = {95, 95, 95}));
    connect(bodyCylinder12.frame_a, revolute7.frame_b) annotation (
      Line(points = {{-146, -182}, {-146, -203.5}, {-174, -203.5}, {-174, -203}}, color = {95, 95, 95}));
    connect(fixedTranslation5.frame_b, revolute8.frame_a) annotation (
      Line(points = {{-300, -222}, {-300, -224}, {-286, -224}}, color = {95, 95, 95}));
    connect(fixedTranslation1.frame_b, fixedTranslation5.frame_a) annotation (
      Line(points = {{-188, 12}, {-320, 12}, {-320, -222}}, color = {95, 95, 95}));
    connect(bodyCylinder13.frame_a, revolute9.frame_b) annotation (
      Line(points = {{-132, 194}, {-132, 185.188}, {-144, 185.188}, {-144, 184.375}, {-164, 184.375}, {-164, 217}}, color = {95, 95, 95}));
    connect(bodyCylinder14.frame_a, revolute9.frame_b) annotation (
      Line(points = {{-130, 230}, {-130, 240.5}, {-164, 240.5}, {-164, 217}}, color = {95, 95, 95}));
    connect(revolute9.frame_a, fixedTranslation9.frame_b) annotation (
      Line(points={{-186,217},{-198,217},{-198,216}},        color = {95, 95, 95}));
    connect(revolute10.frame_b, fixedTranslation9.frame_a) annotation (
      Line(points = {{-236, 218}, {-223, 218}, {-223, 216}, {-218, 216}}));
    connect(fixedTranslation7.frame_b, revolute10.frame_a) annotation (
      Line(points = {{-278, 212}, {-256, 212}, {-256, 218}}, color = {95, 95, 95}));
    connect(revolute8.frame_b, fixedTranslation6.frame_a) annotation (
      Line(points = {{-266, -224}, {-240.5, -224}, {-240.5, -202}, {-236, -202}}));
    connect(revolute7.frame_a, fixedTranslation6.frame_b) annotation (
      Line(points = {{-196, -203}, {-216, -203}, {-216, -202}}, color = {95, 95, 95}));
    connect(revolute6.frame_a, fixedTranslation4.frame_b) annotation (
      Line(points = {{-192, -121}, {-208, -121}, {-208, -124}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_b, fixedTranslation10.frame_a) annotation (
      Line(points = {{-190, -16}, {-190, -32}, {-198, -32}, {-198, -48}}, color = {95, 95, 95}));
    connect(fixedTranslation8.frame_a, revolute5.frame_b) annotation (
      Line(points={{-298,-84},{-298,-114},{-248,-114},{-248,-125}}));
    connect(springDamperParallel.frame_b, fixedTranslation8.frame_b) annotation (
      Line(points = {{-264, -68}, {-252, -68}, {-252, -84}, {-278, -84}}, color = {95, 95, 95}));
    connect(springDamperParallel.frame_a, fixedTranslation10.frame_b) annotation (
      Line(points = {{-264, -48}, {-218, -48}}, color = {95, 95, 95}));
    connect(fixedTranslation12.frame_a, fixedTranslation1.frame_b) annotation (
      Line(points = {{-366, -10}, {-220, -10}, {-220, 12}, {-188, 12}}));
    connect(fixedTranslation12.frame_b, springDamperParallel1.frame_a) annotation (
      Line(points = {{-386, -10}, {-401, -10}, {-401, -72}}, color = {95, 95, 95}));
    connect(fixedTranslation11.frame_a, revolute8.frame_b) annotation (
      Line(points = {{-378, -182}, {-266, -182}, {-266, -224}}, color = {95, 95, 95}));
    connect(fixedTranslation11.frame_b, springDamperParallel1.frame_b) annotation (
      Line(points = {{-398, -182}, {-401, -182}, {-401, -94}}, color = {95, 95, 95}));
    connect(fixedTranslation14.frame_a, fixedTranslation1.frame_b) annotation (
      Line(points = {{-340, 30}, {-188, 30}, {-188, 12}}, color = {95, 95, 95}));
    connect(fixedTranslation14.frame_b, springDamperParallel2.frame_a) annotation (
      Line(points = {{-360, 30}, {-410, 30}, {-410, 46}}, color = {95, 95, 95}));
    connect(revolute10.frame_b, fixedTranslation13.frame_a) annotation (
      Line(points = {{-236, 218}, {-224, 218}, {-224, 152}, {-250, 152}}));
    connect(fixedTranslation13.frame_b, springDamperParallel2.frame_b) annotation (
      Line(points = {{-250, 132}, {-410, 132}, {-410, 66}}));
    connect(bodyCylinder.frame_b, contact_sprocket.frame_b) annotation (
      Line(points={{-110,276},{16.2,276},{16.2,238.4}}));
    connect(speed.flange, revolute.axis) annotation (
      Line(points = {{-450, 268}, {-254, 268}, {-254, 290}, {-186, 290}, {-186, 280}}));
    connect(ramp.y, speed.w_ref) annotation (
      Line(points={{-505.2,268},{-472,268}},    color = {0, 0, 127}));
    connect(revolute9.frame_b, contact.frame_b) annotation (
      Line(points={{-164,217},{-72,217},{-72,175.6},{-35.2,175.6}}));
    connect(contact.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{52,174.8},{202,174.8},{202,65.66}}, color = {95, 95, 95}, thickness = 0.5));
    connect(revolute4.frame_b, contact1.frame_b) annotation (
      Line(points={{-158,127},{-60,127},{-60,79.52},{-38.84,79.52}}, color = {95, 95, 95}));
    connect(contact1.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{44,78.76},{202,78.76},{202,65.66}},
                                                      thickness = 0.5));
    connect(revolute1.frame_b, contact2.frame_b) annotation (
      Line(points={{-156,55},{-58,55},{-58,-8.56},{-34.48,-8.56}}));
    connect(contact2.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{44,-9.28},{202,-9.28},{202,65.66}}, thickness = 0.5));
    connect(revolute2.frame_b, contact3.frame_b) annotation (
      Line(points={{-146,11},{-76,11},{-76,-102.24},{-61.92,-102.24}}, color = {95, 95, 95}));
    connect(contact3.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{34,-103.12},{202,-103.12},{202,65.66}},
                                                          color = {95, 95, 95}, thickness = 0.5));
    connect(revolute3.frame_b, contact4.frame_b) annotation (
      Line(points={{-162,-67},{-82,-67},{-82,-191.44},{-67.02,-191.44}}));
    connect(contact4.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{18,-192.22},{202,-192.22},{202,65.66}},
                                                          thickness = 0.5));
    connect(revolute6.frame_b, contact5.frame_b) annotation (
      Line(points={{-170,-121},{-106,-121},{-106,-286.24},{-69.92,-286.24}}, color = {95, 95, 95}));
    connect(contact5.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{26,-287.12},{202,-287.12},{202,65.66}},
                                                          color = {95, 95, 95}, thickness = 0.5));
    connect(revolute7.frame_b, contact6.frame_b) annotation (
      Line(points={{-174,-203},{-224,-203},{-224,-345.04},{-186.82,-345.04}}));
    connect(contact6.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{-80,-346.02},{202,-346.02},{202,65.66}},
                                                           thickness = 0.5));
    connect(contact_sprocket.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{38,238.2},{124,238.2},{124,65.66},{202,65.66}},thickness = 0.5));
    connect(fixedTranslation15.frame_a, fixedTranslation2.frame_b) annotation (
      Line(points = {{-562, 42}, {-586, 42}, {-586, -16}, {-190, -16}}));
    connect(fixedTranslation15.frame_b, revolute.frame_a) annotation (
      Line(points = {{-562, 62}, {-448, 62}, {-448, 270}, {-196, 270}}, color = {95, 95, 95}));
    connect(revolute4.frame_a, fixedTranslation.frame_a) annotation (
      Line(points={{-180,127},{-242,127},{-242,46},{-208,46}},          color = {95, 95, 95}));
    connect(fixedTranslation.frame_a, fixedTranslation7.frame_a) annotation (
      Line(points = {{-208, 46}, {-320, 46}, {-320, 212}, {-298, 212}}, color = {95, 95, 95}));
    connect(cinghiafinal.frame_a, contact_t.frame_b) annotation (
      Line(points={{336,61.64},{420,61.64},{420,62},{503,62}},
                                            color = {95, 95, 95}, thickness = 0.5));
    connect(force.frame_b, bodyCylinder.frame_a) annotation (
      Line(points = {{-196, 326}, {-130, 326}, {-130, 276}}));
    connect(force.frame_a, world.frame_b) annotation (
      Line(points = {{-216, 326}, {-376, 326}, {-376, 174}, {-496, 174}}));
    connect(realExpression1.y, force.force[1]) annotation (
      Line(points={{-603,332},{-212,332},{-212,337.333}},    color = {0, 0, 127}));
    connect(realExpression.y, force.force[2]) annotation (
      Line(points={{-711,302},{-212,302},{-212,338}},        color = {0, 0, 127}));
    connect(ramp1.y, force.force[3]) annotation (
      Line(points={{-819,166},{-436,166},{-436,372},{-212,372},{-212,338.667}},        color = {0, 0, 127}));
    annotation (
      Icon(coordinateSystem(grid = {2, 0})),
      experiment(StartTime = 0, StopTime = 3, Tolerance = 0.0001, Interval = 0.006));
  end trackterrain;

  model contact_t
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b[80] annotation (
      Placement(visible = true, transformation(origin = {-118, 3}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-118, 4}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    AGILEX_BUNKER_PRO.contact_terrain_new contact_terrain[80] annotation (
        Placement(visible=true, transformation(
          origin={-2,0},
          extent={{-10,-10},{10,10}},
          rotation=0)));
  equation
    for n in 1:80 loop
      connect(contact_terrain[n].frame_b, frame_b[n]);
    end for;
    annotation (
      Icon(coordinateSystem(grid = {2, 0})),
      Diagram(coordinateSystem(grid = {2, 3})));
  end contact_t;

  model cinghiaterrain
    parameter Real start[3] = {0, 0, 0};
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b[80] annotation (
      Placement(visible = true, transformation(origin = {-100, -3}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-100, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a[80] annotation (
      Placement(visible = true, transformation(origin = {100, -9}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {100, -8}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    AGILEX_BUNKER_PRO.element link[80](
      startX={0.36773947,0.3963853,0.42503104,0.4536766,0.4823221,0.51096743,
          0.53961265,0.5682577,0.59690255,0.6255469,0.6541882,0.6828032,
          0.7109103,0.7358329,0.7590959,0.7820636,0.80498064,0.8278859,
          0.8507788,0.87361664,0.8961227,0.9163822,0.92863756,0.92972237,
          0.9196415,0.8999216,0.8736118,0.8453396,0.81690764,0.7884543,
          0.75999576,0.7315221,0.70296574,0.67431283,0.6456575,0.61700207,
          0.58834666,0.55969125,0.53103584,0.50238043,0.47372502,0.44506964,
          0.41641426,0.38775888,0.35910353,0.33044818,0.30179286,0.27313754,
          0.24448223,0.21582696,0.18717167,0.15851632,0.12986113,0.101235166,
          0.07392505,0.048915386,0.024377355,-7.3567404e-05,-0.024447307,-0.0483971,
          -0.069297664,-0.08290987,-0.08756612,-0.08285363,-0.06925454,-0.048120078,
          -0.022466758,0.0042500696,0.031119948,0.058016557,0.08492475,
          0.111871324,0.13903299,0.16720767,0.19585837,0.22450773,0.2531547,
          0.28180113,0.31044737,0.33909348} + fill(start[1], 80),
      startY={-0.14657314,-0.14593814,-0.14529832,-0.1446529,-0.14400178,-0.1433449,
          -0.1426822,-0.14201333,-0.14133595,-0.14063577,-0.13982418,-0.13835242,
          -0.13279773,-0.118662275,-0.1019335,-0.0848015,-0.067601725,-0.050386235,
          -0.033154298,-0.01584938,0.0018847833,0.022145502,0.048039787,
          0.07666748,0.10348362,0.124264404,0.13560154,0.14026824,0.1438363,
          0.14722961,0.15057907,0.15379718,0.15617196,0.15655196,0.15661812,
          0.15663898,0.15665735,0.15668009,0.15670827,0.15674207,0.15678151,
          0.1568266,0.15687735,0.15693374,0.15699577,0.15706347,0.15713681,
          0.15721577,0.15730028,0.15738966,0.15747988,0.15754566,0.15742835,
          0.1561337,0.14748102,0.13349974,0.11870615,0.10376906,0.088706374,
          0.07297853,0.053384025,0.028176269,-9.1155656e-05,-0.028349223,-0.053563982,
          -0.07290421,-0.085662924,-0.09601373,-0.10596052,-0.11583476,-0.12567736,
          -0.13541436,-0.14453399,-0.14972943,-0.15007406,-0.14962608,-0.14904532,
          -0.14843833,-0.14782232,-0.14720018} + fill(start[2], 80),
      startZ=fill(0.0525, 80) + fill(start[3], 80),
      thZ={0.022163361,0.022331439,0.022526653,0.022726301,0.022926996,
          0.023129912,0.023345927,0.023642989,0.024442146,0.028351976,
          0.051533554,0.19409657,0.51524556,0.6233287,0.64086413,0.6438329,
          0.64452124,0.64524186,0.6484496,0.6674835,0.78493977,1.1276712,
          1.5318708,1.929349,2.3289232,2.7331738,2.9777532,3.016711,3.0228894,
          3.02444,3.02908,3.0581698,3.1282597,3.1392727,3.1408634,3.1409514,
          3.1407998,3.1406097,3.1404138,3.1402168,3.1400197,3.1398225,3.1396253,
          3.1394281,3.139231,3.1390338,3.1388373,3.1386445,3.1384742,3.1384454,
          3.1393044,-3.1374555,-3.0961285,-2.8360512,-2.6320481,-2.599088,-2.5931966,
          -2.5880237,-2.5603452,-2.3891828,-2.0668128,-1.7349125,-1.406404,-1.0770421,
          -0.7422006,-0.46210933,-0.36971548,-0.35455695,-0.3518453,-0.35066423,
          -0.34672254,-0.32378414,-0.18343362,-0.012204683,0.015607214,
          0.020264655,0.021184402,0.021500124,0.021714067,0.021885095})
      annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));

  equation
    for n in 1:79 loop
      connect(link[n].frame_b, link[n + 1].frame_a);
    end for;
    connect(link[80].frame_b, link[1].frame_a);
    for n in 1:80 loop
      connect(link[n].frame_b1, frame_b[n]);
      connect(link[n].frame_a1, frame_a[n]);
    end for;
    annotation (
      Icon(coordinateSystem(grid = {2, 0})),
      Diagram(coordinateSystem(grid = {2, 3})),
      Icon(coordinateSystem(grid = {2, 0})),
      Diagram(coordinateSystem(grid = {2, 3})));
  end cinghiaterrain;

  model Surface
    extends Modelica.Mechanics.MultiBody.Visualizers.Advanced.Surface;
    //"implemetation of the surface class based on interoplation between four elevation points with eight paritial derivatives."
    //import SI = Modelica.Units.SI;
    //import MB = Modelica.Mechanics.MultiBody;
    parameter Boolean FlatSurface = true;
    //"If true simpler equations for the functions can be used, false enables uneven surfaces";
    parameter Boolean visSurface = true;
    //"= true if the surfrace shall be shown in the animation window";
    parameter Integer nu1(final min = 2) = 4;
    //"Number of grid points in x direction";
    parameter Integer nv1(final min = 2) = 4;
    //"Number of grid points in z direction";
    parameter Integer IP(final min = 0) = 5;
    //"Number of interpolation points between the grid points (just for visualization)";
    parameter Real PNG = 1;
    // "Filename of the that shall be used as texture";
    //parameter MB.Types.Color Color = {192, 192, 192};
    // "Surface color (mixed with texture)";
    //parameter MB.Types.SpecularCoefficient SpecularCoefficient = 0.1;
    // "Specular coefficient of the road surface without texture";
    parameter Modelica.Units.SI.Length LengthX = 50;
    //"Length of the surface area in x direction";
    parameter Modelica.Units.SI.Length LengthZ = 50;
    //"Length of the surface area in z direction";
    parameter Modelica.Units.SI.Length OffsetX = -LengthX / 2;
    // "Offset in x direction";
    parameter Modelica.Units.SI.Length OffsetZ = -LengthZ / 2;
    // "Offset in z direction";
    parameter Modelica.Units.SI.Position yel[nu1, nv1] = zeros(nu1, nv1);
    //"Matrix with the y (elevation) values at the grid points"
    parameter Real dy_dx[nu1, nv1] = zeros(nu1, nv1);
    //"Matrix with the slope of the grid points in x direction"
    parameter Real dy_dz[nu1, nv1] = zeros(nu1, nv1);
    //"Matrix with the slope of the grid points in y direction"
    parameter Real mu[nu1 - 1, nv1 - 1] = ones(nu1 - 1, nv1 - 1);
    //"Matrix with the friction coefficient for every rectangle";
    final parameter Integer nuIP = nu1 + (nu1 - 1) * IP;
    // "Number of the grid values in x direction with interpolation";
    final parameter Integer nvIP = nv1 + (nv1 - 1) * IP;
    // "Number of the grid values in y direction with interpolation";
    input Modelica.Units.SI.Position d[nu1,nv1];
    function get_eN = get_eN_protected(FlatSurface = FlatSurface, yMatrix = yel, dy_dxMatrix = dy_dx, dy_dzMatrix = dy_dz, nu = nu1, nv = nv1, LengthX = LengthX, LengthZ = LengthZ, OffsetX = OffsetX, OffsetZ = OffsetZ) annotation (
      Documentation(info = "<html>
  
  See the documentation <a href=\"Modelica://WheelsAndTires.Environment.Surface\">Surface</a> for more information.
  
  </html>"));
    function get_elevation = get_elevation_protected(FlatSurface = FlatSurface, yMatrix = yel, dy_dxMatrix = dy_dx, dy_dzMatrix = dy_dz, nu = nu1, nv = nv1, LengthX = LengthX, LengthZ = LengthZ, OffsetX = OffsetX, OffsetZ = OffsetZ) annotation (
      Documentation(info = "<html>
  
  See the documentation <a href=\"Modelica://WheelsAndTires.Environment.Surface\">Surface</a> for more information.
  
  </html>"));
    function get_mu = get_mu_protected(muMatrix = mu, nu = nu1, nv = nv1, LengthX = LengthX, LengthZ = LengthZ, OffsetX = OffsetX, OffsetZ = OffsetZ) annotation (
      Documentation(info = "<html>
  
  See the documentation <a href=\"Modelica://WheelsAndTires.Environment.Surface\">Surface</a> for more information.
  
  </html>"));
    function get_rect = get_rectangle(nu = nu1, nv = nv1, LengthX = LengthX, LengthZ = LengthZ, OffsetX = OffsetX, OffsetZ = OffsetZ);

    function SurfaceFunction
       extends
        Modelica.Mechanics.MultiBody.Interfaces.partialSurfaceCharacteristic;
     input Real x_min "Minimum value of x";
     input Real x_max "Maximum value of x";
     input Real z_min "Minimum value of z";
     input Real z_max "Maximum value of z";
     input Modelica.Units.SI.Distance yMatrix[:, :];
      input Real dy_dxMatrix[:, :];
      input Real dy_dzMatrix[:, :];
     //input Real wz "Factor for angular frequency";
    protected
     Real aux_x;
     //Real A=(z_max-z_min)/2;
    algorithm
      for
      i in 1:nu loop
        aux_x := x_min + (x_max - x_min) * (i - 1) / (nu - 1);
        for j in 1:nv loop
          X[i, j] := aux_x;
          Z[i, j] := z_min + (z_max - z_min) * (j - 1) / (nv - 1);
          Y[i, j] := get_elevation(X[i, j], Z[i, j], yMatrix = yMatrix, dy_dxMatrix = dy_dxMatrix, dy_dzMatrix = dy_dzMatrix);
        end for;
      end for;
      //y[i,j];
      if multiColoredSurface then
        C := {{((-Y[i, j]) + 0) * 200000, 255, 0} for j in 1:nv, i in 1:nu};
      end if;
    end SurfaceFunction;

   redeclare function surfaceCharacteristic = SurfaceFunction(x_min=OffsetX + 1e-6, x_max=LengthX+OffsetX - 1e-6, z_min=OffsetZ + 1e-6, z_max=LengthZ+OffsetZ - 1e-6, yMatrix = yel+d, dy_dxMatrix = dy_dx+d, dy_dzMatrix = dy_dz+d);

    function get_rectangle
      //"finds the rectangle which calculates the elevation and normal vector for the current position of the tire."
      input Modelica.Units.SI.Position x;
      //"acutal position to find active recthangle and relative position in x direction";
      input Modelica.Units.SI.Position z;
      //"acutal position to find active recthangle and relative position in z direction";
      input Integer nu;
      // "Number of grid points in x direction";
      input Integer nv;
      // "Number of grid points in z direction";
      input Modelica.Units.SI.Length LengthX;
      // "Length of the surface area in x direction";
      input Modelica.Units.SI.Length LengthZ;
      // "Length of the surface area in z direction";
      input Modelica.Units.SI.Distance OffsetX;
      // "Offset of the surface area in x direction";
      input Modelica.Units.SI.Distance OffsetZ;
      // "Offset of the surface area in z direction";
      output Integer rectX;
      // "index of the active rectange in x direction";
      output Integer rectZ;
      // "index of the active rectange in z direction";
      output Boolean inArea;
      //"true: tire is in the defined surface area, false: it is not.";
      output Modelica.Units.SI.Distance localX;
      //"relative position in the active rectangle in x direction.";
      output Modelica.Units.SI.Distance localZ;
      //"relative position in the active rectangle in z direction.";
    protected
      Modelica.Units.SI.Distance rectangleX;
      Modelica.Units.SI.Distance rectangleZ;
    algorithm
      rectX := integer(ceil((x - OffsetX) / LengthX * (nu - 1)));
      rectZ := integer(ceil((z - OffsetZ) / LengthZ * (nv - 1)));
      if rectX > nu - 1 or rectX < 1 or rectZ > nv - 1 or rectZ < 1 then
        inArea := false;
      else
        inArea := true;
      end if;
      rectangleX := LengthX / (nu - 1);
      rectangleZ := LengthZ / (nv - 1);
      // getting local position (in the rectangle the tire is in at the time of calling)
      localX := mod(x - LengthX - OffsetX, rectangleX);
      localZ := mod(z - LengthZ - OffsetZ, rectangleZ);
    end get_rectangle;

  protected
    function get_eN_protected
      input Modelica.Units.SI.Position x;
      //"actual position to get the normal vector eN in x direction."
      input Modelica.Units.SI.Position z;
      //"actual position to get the normal vector eN in z direction."
      input Boolean FlatSurface;
      input Modelica.Units.SI.Distance yMatrix[:, :];
      input Real dy_dxMatrix[:, :];
      input Real dy_dzMatrix[:, :];
      input Integer nu;
      input Integer nv;
      input Modelica.Units.SI.Length LengthX;
      input Modelica.Units.SI.Length LengthZ;
      input Modelica.Units.SI.Distance OffsetX;
      input Modelica.Units.SI.Distance OffsetZ;
      output Real eN[3];
      // "normal vector on the surface.";
    protected
      Integer xIndex;
      Integer zIndex;
      Modelica.Units.SI.Distance localX;
      Modelica.Units.SI.Distance localZ;
      Boolean inArea;
      Modelica.Units.SI.Distance rectangleY[2, 2];
      Real rectangleDy_dx[2, 2];
      Real rectangleDy_dz[2, 2];
      Real X[4];
      Real Z[4];
      Real G[4, 4];
      constant Real H[4, 4] = [2, -3, 0, 1; -2, 3, 0, 0; 1, -2, 1, 0; 1, -1, 0, 0];
      constant Real dH[4, 4] = [0, 6, -6, 0; 0, -6, 6, 0; 0, 3, -4, 1; 0, 3, -2, 0];
      Real dy_dx;
      Real dy_dz;
      Modelica.Units.SI.Distance factorX;
      Modelica.Units.SI.Distance factorZ;
      Modelica.Units.SI.Length localXNorm;
      Modelica.Units.SI.Length localZNorm;
      Modelica.Units.SI.Length rectangleX;
      Modelica.Units.SI.Length rectangleZ;
      Real dN[3];
    algorithm
      (xIndex, zIndex, inArea, localX, localZ) := get_rectangle(x = x, z = z, nu = nu, nv = nv, LengthX = LengthX, LengthZ = LengthZ, OffsetX = OffsetX, OffsetZ = OffsetZ);
      if inArea and not FlatSurface then
        rectangleX := LengthX / (nu - 1);
        rectangleZ := LengthZ / (nv - 1);
        rectangleY := yMatrix[xIndex:xIndex + 1, zIndex:zIndex + 1];
        rectangleDy_dx := dy_dxMatrix[xIndex:xIndex + 1, zIndex:zIndex + 1];
        rectangleDy_dz := dy_dzMatrix[xIndex:xIndex + 1, zIndex:zIndex + 1];
        G := [rectangleY, rectangleDy_dz * rectangleZ; rectangleDy_dx * rectangleX, zeros(2, 2)];
        localXNorm := localX / rectangleX;
        localZNorm := localZ / rectangleZ;
        X := {localXNorm ^ 3, localXNorm ^ 2, localXNorm, 1};
        Z := {localZNorm ^ 3, localZNorm ^ 2, localZNorm, 1};
        dy_dx := dH * X * G * (H * Z) / rectangleX;
        dy_dz := H * X * G * (dH * Z) / rectangleZ;
        dN := cross({0, dy_dz, 1}, {1, dy_dx, 0});
        eN := dN / sqrt(dN * dN);
      else
        eN := {0, 1, 0};
      end if;
    end get_eN_protected;

    function get_elevation_protected
      input Modelica.Units.SI.Position x;
      //"actual position to get the elevation of the surface (y-coordinate) in x direction"
      input Modelica.Units.SI.Position z;
      //"actual position to get the elevation of the surface (y-coordinate) in z direction"
      input Boolean FlatSurface;
      input Modelica.Units.SI.Distance yMatrix[:, :];
      input Real dy_dxMatrix[:, :];
      input Real dy_dzMatrix[:, :];
      input Integer nu;
      input Integer nv;
      input Modelica.Units.SI.Length LengthX;
      input Modelica.Units.SI.Length LengthZ;
      input Modelica.Units.SI.Distance OffsetX;
      input Modelica.Units.SI.Distance OffsetZ;
      output Modelica.Units.SI.Position elevation;
    protected
      Integer xIndex;
      Integer zIndex;
      Modelica.Units.SI.Distance localX;
      Modelica.Units.SI.Distance localZ;
      Boolean inArea;
      Modelica.Units.SI.Distance rectangleY[2, 2];
      Real rectangleDy_dx[2, 2];
      Real rectangleDy_dz[2, 2];
      Real X[4];
      Real Z[4];
      Real G[4, 4];
      Modelica.Units.SI.Distance rectangleX;
      Modelica.Units.SI.Distance rectangleZ;
      Real localXNorm;
      Real localZNorm;
      constant Real H[4, 4] = [2, -3, 0, 1; -2, 3, 0, 0; 1, -2, 1, 0; 1, -1, 0, 0];
      //"Matrix used for the interpolation";
    algorithm
      (xIndex, zIndex, inArea, localX, localZ) := get_rectangle(x = x, z = z, nu = nu, nv = nv, LengthX = LengthX, LengthZ = LengthZ, OffsetX = OffsetX, OffsetZ = OffsetZ);
      if inArea and not FlatSurface then
        rectangleX := LengthX / (nu - 1);
        rectangleZ := LengthZ / (nv - 1);
        rectangleY := yMatrix[xIndex:xIndex + 1, zIndex:zIndex + 1];
        rectangleDy_dx := dy_dxMatrix[xIndex:xIndex + 1, zIndex:zIndex + 1];
        rectangleDy_dz := dy_dzMatrix[xIndex:xIndex + 1, zIndex:zIndex + 1];
        G := [rectangleY, rectangleDy_dz * rectangleZ; rectangleDy_dx * rectangleX, zeros(2, 2)];
        localXNorm := localX / rectangleX;
        localZNorm := localZ / rectangleZ;
        X := {localXNorm ^ 3, localXNorm ^ 2, localXNorm, 1};
        Z := {localZNorm ^ 3, localZNorm ^ 2, localZNorm, 1};
        elevation := H * X * G * (H * Z);
      else
        elevation := 0;
      end if;
    end get_elevation_protected;

    function get_mu_protected
      input Modelica.Units.SI.Position x;
      //"actual position to get the friction coefficient of the surface in x direction"
      input Modelica.Units.SI.Position z;
      //"actual position to get the friction coefficient of the surface in z direction"
      input Real muMatrix[:, :];
      //additional input(s)
      input Integer nu;
      input Integer nv;
      input Modelica.Units.SI.Length LengthX;
      input Modelica.Units.SI.Length LengthZ;
      input Modelica.Units.SI.Distance OffsetX;
      input Modelica.Units.SI.Distance OffsetZ;
      output Real mu;
      //"the friction coefficient";
    protected
      Integer xIndex;
      Integer zIndex;
      Boolean inArea;
    algorithm
      (xIndex, zIndex, inArea) := get_rectangle(x = x, z = z, nu = nu, nv = nv, LengthX = LengthX, LengthZ = LengthZ, OffsetX = OffsetX, OffsetZ = OffsetZ);
      if inArea then
        mu := muMatrix[xIndex, zIndex];
      else
        mu := 1;
      end if;
    end get_mu_protected;

  end Surface;

  model contact_terrain_new
    parameter Integer nx = 1;
    parameter Integer nz = 1;
    parameter Modelica.Units.SI.Length len = 2.28/80;
    parameter Modelica.Units.SI.Length lar = 0.15;
    parameter Modelica.Units.SI.Length thick = 0.05;
    parameter Real k = 1646666.667;
    //290000;
    parameter Real kv = 5e4;
    //290000;
    parameter Real n = 0.13;
    ///2;
    parameter Real c = 70000;
    //70000 * 0.5;
    parameter Real tanf = 0.67;
    parameter Real K = 0.02;
    parameter Real gammas = 15000;
    parameter Real Nf = tan(Modelica.Constants.pi / 4 + 0.5 * 0.590307) ^ 2;
    // parameter Real w = 10;
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b annotation (
      Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-100, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Units.SI.Position dr[nx * nz, 3],dr1[3],dr2[3];
    Modelica.Units.SI.Position d[nx * nz, 3];
    Modelica.Units.SI.Velocity vjx[nx * nz], vjz[nx * nz], modv[nx * nz];
    //(each stateSelect=StateSelect.prefer)
    Modelica.Units.SI.Position jx[nx * nz], jz[nx * nz], j[nx * nz];
  //, modj[nx * nz];
    Real cos[nx * nz], sin[nx * nz];
    //, cos1[nx * nz], sin1[nx * nz];
    Modelica.Units.SI.Distance deltaZ[nx * nz];
    Modelica.Units.SI.Pressure p[nx * nz];
    Modelica.Units.SI.Pressure t[nx * nz];
    Modelica.Units.SI.Force Fbr[nx];
    Modelica.Units.SI.Force Fbl[nx];
    //, cos1[nx * nz], sin1[nx * nz];
    //Real smooth11[nx], smooth12[nx];
    Modelica.Units.SI.Force Fn, Ft, Fl;
    Modelica.Units.SI.Torque T1, T2, T3;
    Modelica.Units.SI.Position dd[Surface.nu1,Surface.nv1];
    Integer xIndex,zIndex, G[4,4],H[4,4];
    Boolean inArea;
    Modelica.Units.SI.Distance localX;
    Modelica.Units.SI.Distance localZ;
    Real X[4];
    Real Z[4];
    Modelica.Units.SI.Distance rectangleX;
    Modelica.Units.SI.Distance rectangleZ;
    Real localXNorm;
    Real localZNorm;
    Boolean inArea1;
    Modelica.Units.SI.Distance localX1;
    Modelica.Units.SI.Distance localZ1;
    Real X1[4];
    Real Z1[4];
    Real localXNorm1;
    Real localZNorm1;
    Integer xIndex1,zIndex1;
    Boolean inArea2;
    Modelica.Units.SI.Distance localX2;
    Modelica.Units.SI.Distance localZ2;
    Real X2[4];
    Real Z2[4];
    Real localXNorm2;
    Real localZNorm2;
    Integer xIndex2,zIndex2;
    function smooth1
      input Modelica.Units.SI.Distance d1;
      input Real w = 1e4;
      output Real s;
    algorithm
      //if d1 == 0 then
      //s := 0;
      //else
      s := (Modelica.Math.tanh(w * (d1 - 1e-3)) + 1) / 2;
      //end if;//(Modelica.Math.tanh(w * (d1 - 1e-3)) + 1) / 2;//((Modelica.Math.atan(w * (d1 - 1e-3)) + Modelica.Constants.pi / 2) / Modelica.Constants.pi);//
      //((Modelica.Math.atan(w * d1) + Modelica.Constants.pi / 2) / Modelica.Constants.pi);
    end smooth1;

    outer AGILEX_BUNKER_PRO.Surface Surface;
  equation
  for n in 1:nx loop
      for m in 1:nz loop
        dr[(n - 1) * nz + m, :] = frame_b.r_0 + Modelica.Mechanics.MultiBody.Frames.resolve1(frame_b.R, {len * (n / nx) - len * (0.5 / nx), 0, lar * (m / nz) - lar * (0.5 / nz)});
        d[(n - 1) * nz + m, :] = {len * (n / nx) - len * (0.5 / nx), 0, lar * (m / nz) - lar * (0.5 / nz)};
        if noEvent(dr[(n - 1) * nz + m, 2] - Surface.get_elevation(dr[(n - 1) * nz + m, 1], dr[(n - 1) * nz + m, 3]) < 0) then
          vjx[(n - 1) * nz + m] = der(dr[(n - 1) * nz + m, :]) * Modelica.Mechanics.MultiBody.Frames.resolve1(frame_b.R, {1, 0, 0});
          vjz[(n - 1) * nz + m] = der(dr[(n - 1) * nz + m, :]) * Modelica.Mechanics.MultiBody.Frames.resolve1(frame_b.R, {0, 0, 1});
          der(jx[(n - 1) * nz + m]) = vjx[(n - 1) * nz + m];
          der(jz[(n - 1) * nz + m]) = vjz[(n - 1) * nz + m];
          deltaZ[(n - 1) * nz + m] = -(dr[(n - 1) * nz + m, 2] - Surface.get_elevation(dr[(n - 1) * nz + m, 1], dr[(n - 1) * nz + m, 3]));
          p[(n - 1) * nz + m] = -smooth(0, max(0, k * deltaZ[(n - 1) * nz + m] ^ n - kv * der(dr[(n - 1) * nz + m, 2])));
        else
          vjx[(n - 1) * nz + m] = 0;
          vjz[(n - 1) * nz + m] = 0;
          der(jx[(n - 1) * nz + m]) = -1000 * jx[(n - 1) * nz + m];
          der(jz[(n - 1) * nz + m]) = -1000 * jz[(n - 1) * nz + m];
          deltaZ[(n - 1) * nz + m] = 0;
          p[(n - 1) * nz + m] = 0;
        end if;
        j[(n - 1) * nz + m] = sqrt(jx[(n - 1) * nz + m] ^ 2 + jz[(n - 1) * nz + m] ^ 2);
        //sin1[(n - 1) * nz + m] = vjz[(n - 1) * nz + m] / (Modelica.Constants.small + sqrt(vjx[(n - 1) * nz + m] ^ 2 + vjz[(n - 1) * nz + m] ^ 2));
        //cos1[(n - 1) * nz + m] = vjx[(n - 1) * nz + m] / (Modelica.Constants.small + sqrt(vjx[(n - 1) * nz + m] ^ 2 + vjz[(n - 1) * nz + m] ^ 2));
        modv[(n - 1) * nz + m] = smooth(0, max(1e-3, sqrt(vjx[(n - 1) * nz + m] ^ 2 + vjz[(n - 1) * nz + m] ^ 2)));
        sin[(n - 1) * nz + m] = vjz[(n - 1) * nz + m] / modv[(n - 1) * nz + m];
        cos[(n - 1) * nz + m] = vjx[(n - 1) * nz + m] / modv[(n - 1) * nz + m];
        //p[(n - 1) * nz + m] = -max(0, k * deltaZ[(n - 1) * nz + m] ^ n + kv * der(deltaZ[(n - 1) * nz + m]));
        t[(n - 1) * nz + m] = (c - tanf * p[(n - 1) * nz + m]) * (1 - exp(-j[(n - 1) * nz + m] / K));
        //modj[(n - 1) * nz + m] = max(1e-6, j[(n - 1) * nz + m]);
        //sin[(n - 1) * nz + m] = jz[(n - 1) * nz + m] / modj[(n - 1) * nz + m];
        //cos[(n - 1) * nz + m] = jx[(n - 1) * nz + m] / modj[(n - 1) * nz + m];
        //cos[(n - 1) * nz + m] = (Modelica.Math.tanh(1e3 * vjx[(n - 1) * nz + m])) * sqrt(1 - sin[(n - 1) * nz + m]^2);
      end for;
      Fbr[n] = (0.5 * gammas * deltaZ[nz * n] ^ 2 * Nf + 2 * c * deltaZ[nz * n] * sqrt(Nf)) * (len / nx) * smooth1(vjz[nz * n]);
      Fbl[n] = (0.5 * gammas * deltaZ[nz * (n - 1) + 1] ^ 2 * Nf + 2 * c * deltaZ[nz * (n - 1) + 1] * sqrt(Nf)) * (len / nx) * (-smooth1(-vjz[nz * (n - 1) + 1]));
    end for;
    Ft = t * cos * (len * lar / (nx * nz));
    Fn = sum(p) * (len * lar / (nx * nz));
    Fl = t * sin * (len * lar / (nx * nz)) + sum(Fbr) + sum(Fbl);
    frame_b.f = {Ft, Fn, Fl};
    //frame_b.f = {t * cos * (len * lar / (nx * nz)), sum(p) * (len * lar / (nx * nz)), t * sin * (len * lar / (nx * nz)) + sum(Fbr) + sum(Fbl)};
    T1 = -p * d[:, 3] * (len * lar / (nx * nz));
    T2 = sum(t[i] * (len * lar / (nx * nz)) * (d[i, 3] * cos[i] - d[i, 1] * sin[i]) for i in 1:nx * nz) + sum(Fbr[i] * d[nz * i, 1] for i in 1:nx) + sum(Fbl[i] * d[nz * (i - 1) + 1, 1] for i in 1:nx);
    T3 = p * d[:, 1] * (len * lar / (nx * nz));
    frame_b.t = {T1, T2, T3};
    //frame_b.t = {-p * d[:, 3] * (len * lar / (nx * nz)), sum(t[i] * (len * lar / (nx * nz)) * (d[i, 3] * cos[i] - d[i, 1] * sin[i]) for i in 1:nx * nz) + sum(Fbr[i] * d[nz * i, 1] for i in 1:nx) + sum(Fbl[i] * d[nz * (i - 1) + 1, 1] for i in 1:nx), p * d[:, 1] * (len * lar / (nx * nz))};
    //for i in 1:nx loop
    //smooth11[i] = smooth1(vjz[nz * i]);
    //smooth12[i] = -smooth1(-vjz[nz * (i - 1) + 1]);
    //end for;
  algorithm
    dr1[:] :=frame_b.r_0 + Modelica.Mechanics.MultiBody.Frames.resolve1(frame_b.R,
      {len*0.5,0,lar*0.25});
      dr2[:] :=frame_b.r_0 + Modelica.Mechanics.MultiBody.Frames.resolve1(frame_b.R,
      {len*0.5,0,lar*0.75});
    G:=[1,1,1,1;1,1,1,1;1,1,0,0;1,1,0,0];
    H:=[2, -3, 0, 1; -2, 3, 0, 0; 1, -2, 1, 0; 1, -1, 0, 0];
    (xIndex, zIndex, inArea, localX, localZ) := Surface.get_rect(dr1[1],dr1[3]);
    (xIndex1, zIndex1, inArea1, localX1, localZ1) := Surface.get_rect(dr2[1],dr2[3]);

    rectangleX := Surface.LengthX / (Surface.nu1 - 1);
    rectangleZ := Surface.LengthZ / (Surface.nv1 - 1);
    if xIndex==xIndex1 and zIndex==zIndex1 then
    (xIndex2, zIndex2, inArea2, localX2, localZ2) := Surface.get_rect(frame_b.r_0[1],frame_b.r_0[3]);
    localXNorm2 := localX2 / rectangleX;
    localZNorm2 := localZ2 / rectangleZ;
    X2 := {localXNorm2 ^ 3, localXNorm2 ^ 2, localXNorm2, 1};
    Z2 := {localZNorm2 ^ 3, localZNorm2 ^ 2, localZNorm2, 1};
  dd:=fill(0,Surface.nu1,Surface.nv1);
    for n in 1:Surface.nu1 loop
      for m in 1:Surface.nv1 loop
        if (n == xIndex or n == xIndex + 1) and  (m == zIndex or m == zIndex + 1) then
          dd[n,m]:=-(sum(deltaZ)/(nx*nz))/(H * X2 * G * (H * Z2));
        else
          dd[n,m]:=0;
        end if;
        end for;
      end for;
    else
    localXNorm := localX / rectangleX;
    localZNorm := localZ / rectangleZ;
    X := {localXNorm ^ 3, localXNorm ^ 2, localXNorm, 1};
    Z := {localZNorm ^ 3, localZNorm ^ 2, localZNorm, 1};

    localXNorm1 := localX1 / rectangleX;
    localZNorm1 := localZ1 / rectangleZ;
    X1 := {localXNorm1 ^ 3, localXNorm1 ^ 2, localXNorm1, 1};
    Z1 := {localZNorm1 ^ 3, localZNorm1 ^ 2, localZNorm1, 1};
  dd:=fill(0,Surface.nu1,Surface.nv1);
    for n in 1:Surface.nu1 loop
      for m in 1:Surface.nv1 loop
        if (n == xIndex or n == xIndex + 1) and  (m == zIndex or m == zIndex + 1) or  (n == xIndex1 or n == xIndex1 + 1) and  (m == zIndex1 or m == zIndex1 + 1) then
          if (n == xIndex or n == xIndex + 1) and  (m == zIndex or m == zIndex + 1) then
          dd[n,m]:=dd[n,m]-(sum(deltaZ)/(nx*nz))/(H * X * G * (H * Z))*0.5;
          else
          dd[n,m]:=dd[n,m]-(sum(deltaZ)/(nx*nz))/(H * X1 * G * (H * Z1))*0.5;
          end if;
        else
          dd[n,m]:=0;
        end if;
        end for;
      end for;
    end if;

    annotation (
      Icon(coordinateSystem(grid = {2, 0})),
      Diagram(coordinateSystem(grid = {2, 3})));
  end contact_terrain_new;

  model trackterrainfinal
    parameter Real start[3] = {0, 0, 0};
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder(diameter = 0.16,
      length=0.105,
      lengthDirection={0,0,1},                                                                                               r = {0, 0, 0.0525}, r_0(each fixed = true, start = start), useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-118, 278}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute(animation = false, useAxisFlange = true) annotation (
      Placement(visible = true, transformation(origin = {-186, 270}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder4(diameter = 0.088, length = 0.075, r = {0, 0, 0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-124, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute1(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-167, 55}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder1(diameter = 0.088, length = 0.075, r = {0, 0, -0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-124, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder2(diameter = 0.088, length = 0.075, r = {0, 0, 0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-126, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute2(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-155, 13}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder5(diameter = 0.088, length = 0.075, r = {0, 0, -0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-129, -35}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder6(diameter = 0.088, length = 0.075, r = {0, 0, 0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-132, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute3(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-173, -67}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder7(diameter = 0.130, length = 0.075, r = {0, 0, 0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-124, 114}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder8(diameter = 0.130, length = 0.075, r = {0, 0, -0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-124, 160}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute4(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-169, 127}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation(animation = true, r = {-0.153, 0.040, 0}) annotation (
      Placement(visible = true, transformation(origin = {-198, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation1(r = {-0.350, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-188, 22}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation2(r = {-0.25, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-190, -6}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute5(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-261, -125}, extent = {{-13, 13}, {13, -13}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation3(animation = false, r = {0.250, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-300, -126}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation4(r = {-0.192, -0.192, 0}) annotation (
      Placement(visible = true, transformation(origin = {-218, -124}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder9(diameter = 0.105, length = 0.075, r = {0, 0, -0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-132, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute6(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-181, -121}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder10(diameter = 0.105, length = 0.075, r = {0, 0, 0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-128, -148}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder3(diameter = 0.088, length = 0.075, r = {0, 0, -0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-122, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder11(diameter = 0.105, length = 0.075, r = {0, 0, 0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-134, -222}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute7(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-185, -203}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation5(animation = false, r = {0.281, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-310, -222}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder12(diameter = 0.105, length = 0.075, r = {0, 0, -0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-134, -182}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute8(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-276, -224}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation7(r = {-0.087, -0.029, 0}) annotation (
      Placement(visible = true, transformation(origin = {-286, 212}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder13(diameter = 0.105, length = 0.075, r = {0, 0, 0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-122, 194}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute9(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-175, 217}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder14(diameter = 0.105, length = 0.075, r = {0, 0, -0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-118, 230}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation9(r = {-0.071, -0.120, 0}) annotation (
      Placement(visible = true, transformation(origin = {-208, 216}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute10(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-264, 240}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation6(r = {-0.192, -0.192, 0}) annotation (
      Placement(visible = true, transformation(origin = {-226, -202}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation8(animation = false, r = {-0.134, -0.134, 0}) annotation (
      Placement(visible = true, transformation(origin = {-288, -84}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation10(animation = false, r = {0.064, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-208, -48}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel springDamperParallel(
      c=100000,
      d=1000,                                                                                             fixedRotationAtFrame_a = false, fixedRotationAtFrame_b = false,
      s_unstretched=0.15)                                                                                                                                                                       annotation (
      Placement(visible = true, transformation(origin = {-264, -58}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation11(animation = false, r = {-0.134, -0.134, 0}) annotation (
      Placement(visible = true, transformation(origin = {-388, -182}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation12(animation = false, r = {0.098, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-376, -10}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel springDamperParallel1(
      c=100000,
      d=1000,                                                                                              fixedRotationAtFrame_a = false, fixedRotationAtFrame_b = false,
      s_unstretched=0.14)                                                                                                                                                                        annotation (
      Placement(visible = true, transformation(origin = {-401, -83}, extent = {{-11, -11}, {11, 11}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel springDamperParallel2(
      c=100000,
      d=1000,                                                                                              fixedRotationAtFrame_a = false, fixedRotationAtFrame_b = false,
      s_unstretched=0.17)                                                                                                                                                                        annotation (
      Placement(visible = true, transformation(origin = {-410, 56}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation13(animation = false, r = {-0.044, -0.065, 0}) annotation (
      Placement(visible = true, transformation(origin = {-250, 142}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation14(animation = false, r = {0.278, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-350, 30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    AGILEX_BUNKER_PRO.contact_sprocket contact_sprocket annotation (Placement(
          visible=true, transformation(
          origin={28,238},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    AGILEX_BUNKER_PRO.cinghiaterrain cinghiafinal(start=start) annotation (
        Placement(visible=true, transformation(
          origin={271,67},
          extent={{-67,-67},{67,67}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact(R=0.0525) annotation (Placement(visible=
            true, transformation(
          origin={12,174},
          extent={{-40,-40},{40,40}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact1(R=0.065) annotation (Placement(visible=
            true, transformation(
          origin={6,78},
          extent={{-38,-38},{38,38}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact2(R=0.044) annotation (Placement(visible=
            true, transformation(
          origin={8,-10},
          extent={{-36,-36},{36,36}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact3(R=0.044) annotation (Placement(visible=
            true, transformation(
          origin={-10,-104},
          extent={{-44,-44},{44,44}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact4(R=0.044) annotation (Placement(visible=
            true, transformation(
          origin={-21,-193},
          extent={{-39,-39},{39,39}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact5(R=0.0525) annotation (Placement(visible=
            true, transformation(
          origin={-18,-288},
          extent={{-44,-44},{44,44}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact6(R=0.0525) annotation (Placement(visible=
            true, transformation(
          origin={-129,-347},
          extent={{-49,-49},{49,49}},
          rotation=0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation15(r={-0.105,
          -0.105,-0.0525})                                                                              annotation (
      Placement(visible = true, transformation(origin = {-562, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation16(animation = false, r={0.458,
          0.065,0.0525})                                                                                                annotation (
      Placement(visible = true, transformation(origin = {-528, -44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a annotation (
      Placement(visible = true, transformation(origin = {-664, -36}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-100, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b[80] annotation (
      Placement(visible = true, transformation(origin = {494, 64}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {98, 0}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Components.Damper damper(d=7)
      annotation (Placement(transformation(extent={{-198,306},{-178,326}})));
    Modelica.Mechanics.MultiBody.Joints.Prismatic prismatic(useAxisFlange=true)
                       annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-204,88})));
    Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_a annotation (
        Placement(transformation(extent={{-6,90},{14,110}}), iconTransformation(
            extent={{-6,90},{14,110}})));
    Modelica.Mechanics.MultiBody.Sensors.AbsoluteAngularVelocity
      absoluteAngularVelocity1
      annotation (Placement(transformation(extent={{-128,324},{-108,344}})));
    Modelica.Blocks.Continuous.Filter filter(f_cut=1)
      annotation (Placement(transformation(extent={{-86,324},{-66,344}})));
    Modelica.Mechanics.MultiBody.Sensors.AbsoluteAngularVelocity
      absoluteAngularVelocity2
      annotation (Placement(transformation(extent={{0,298},{20,318}})));
    Modelica.Blocks.Continuous.Filter filter1(f_cut=1)
      annotation (Placement(transformation(extent={{42,298},{62,318}})));
    Modelica.Mechanics.MultiBody.Sensors.AbsoluteAngularVelocity
      absoluteAngularVelocity3
      annotation (Placement(transformation(extent={{36,258},{56,278}})));
    Modelica.Blocks.Continuous.Filter filter2(f_cut=1)
      annotation (Placement(transformation(extent={{78,258},{98,278}})));
    Modelica.Mechanics.MultiBody.Sensors.AbsoluteAngularVelocity
      absoluteAngularVelocity4
      annotation (Placement(transformation(extent={{156,192},{176,212}})));
    Modelica.Blocks.Continuous.Filter filter3(f_cut=1)
      annotation (Placement(transformation(extent={{198,192},{218,212}})));
    Modelica.Mechanics.MultiBody.Sensors.AbsoluteAngularVelocity
      absoluteAngularVelocity5
      annotation (Placement(transformation(extent={{316,-66},{336,-46}})));
    Modelica.Blocks.Continuous.Filter filter4(f_cut=1)
      annotation (Placement(transformation(extent={{358,-66},{378,-46}})));
    Modelica.Mechanics.MultiBody.Sensors.AbsoluteAngularVelocity
      absoluteAngularVelocity6
      annotation (Placement(transformation(extent={{254,-92},{274,-72}})));
    Modelica.Blocks.Continuous.Filter filter5(f_cut=1)
      annotation (Placement(transformation(extent={{296,-92},{316,-72}})));
    Modelica.Mechanics.MultiBody.Sensors.AbsoluteAngularVelocity
      absoluteAngularVelocity7
      annotation (Placement(transformation(extent={{254,-182},{274,-162}})));
    Modelica.Blocks.Continuous.Filter filter6(f_cut=1)
      annotation (Placement(transformation(extent={{296,-182},{316,-162}})));
    Modelica.Mechanics.MultiBody.Sensors.AbsoluteAngularVelocity
      absoluteAngularVelocity8
      annotation (Placement(transformation(extent={{266,-246},{286,-226}})));
    Modelica.Blocks.Continuous.Filter filter7(f_cut=1)
      annotation (Placement(transformation(extent={{308,-246},{328,-226}})));
    Modelica.Mechanics.Translational.Sources.Position position(exact=true)
      annotation (Placement(transformation(extent={{-252,98},{-232,118}})));
    Modelica.Blocks.Sources.RealExpression realExpression
      annotation (Placement(transformation(extent={{-278,72},{-258,92}})));
  equation
    connect(revolute.frame_b, bodyCylinder.frame_a) annotation (
      Line(points = {{-176, 270}, {-140, 270}, {-140, 278}, {-128, 278}}));
    connect(bodyCylinder4.frame_a, revolute1.frame_b) annotation (
      Line(points = {{-134, 60}, {-134, 52.75}, {-138, 52.75}, {-138, 55.5}, {-156, 55.5}, {-156, 55}}, color = {95, 95, 95}));
    connect(bodyCylinder1.frame_a, revolute2.frame_b) annotation (
      Line(points = {{-134, 28}, {-146, 28}, {-146, 13}}, color = {95, 95, 95}));
    connect(bodyCylinder2.frame_a, revolute2.frame_b) annotation (
      Line(points = {{-136, -4}, {-136, -7.5}, {-146, -7.5}, {-146, 13}}, color = {95, 95, 95}));
    connect(bodyCylinder6.frame_a, revolute3.frame_b) annotation (
      Line(points = {{-142, -70}, {-142, -75.5}, {-162, -75.5}, {-162, -67}}, color = {95, 95, 95}));
    connect(bodyCylinder5.frame_a, revolute3.frame_b) annotation (
      Line(points = {{-140, -35}, {-140, -29.5}, {-166, -29.5}, {-166, -55.25}, {-162, -55.25}, {-162, -67}}, color = {95, 95, 95}));
    connect(bodyCylinder8.frame_a, revolute4.frame_b) annotation (
      Line(points = {{-134, 160}, {-134, 157.5}, {-158, 157.5}, {-158, 127}}, color = {95, 95, 95}));
    connect(bodyCylinder7.frame_a, revolute4.frame_b) annotation (
      Line(points = {{-134, 114}, {-134, 110.375}, {-158, 110.375}, {-158, 127}}, color = {95, 95, 95}));
    connect(fixedTranslation.frame_b, revolute1.frame_a) annotation (
      Line(points = {{-188, 46}, {-188, 55}, {-178, 55}}));
    connect(fixedTranslation.frame_b, fixedTranslation1.frame_a) annotation (
      Line(points = {{-188, 46}, {-188, 32}}));
    connect(fixedTranslation1.frame_b, revolute2.frame_a) annotation (
      Line(points = {{-188, 12}, {-188, 13}, {-164, 13}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_b, revolute3.frame_a) annotation (
      Line(points = {{-190, -16}, {-191.25, -16}, {-191.25, -32}, {-192.5, -32}, {-192.5, -40}, {-184, -40}, {-184, -67}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_a, fixedTranslation1.frame_b) annotation (
      Line(points = {{-190, 4}, {-189, 4}, {-189, 12}, {-188, 12}}));
    connect(fixedTranslation3.frame_b, revolute5.frame_a) annotation (
      Line(points = {{-290, -126}, {-303, -126}, {-303, -125}, {-274, -125}}, color = {95, 95, 95}));
    connect(revolute5.frame_b, fixedTranslation4.frame_a) annotation (
      Line(points = {{-248, -125}, {-248, -137.5}, {-228, -137.5}, {-228, -124}}));
    connect(bodyCylinder10.frame_a, revolute6.frame_b) annotation (
      Line(points = {{-138, -148}, {-138, -151.813}, {-146, -151.813}, {-146, -149.625}, {-170, -149.625}, {-170, -121}}, color = {95, 95, 95}));
    connect(bodyCylinder9.frame_a, revolute6.frame_b) annotation (
      Line(points = {{-142, -110}, {-142, -106.75}, {-160, -106.75}, {-160, -105.5}, {-170, -105.5}, {-170, -121}}, color = {95, 95, 95}));
    connect(bodyCylinder3.frame_a, revolute1.frame_b) annotation (
      Line(points = {{-132, 86}, {-156, 86}, {-156, 55}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_b, fixedTranslation3.frame_a) annotation (
      Line(points = {{-190, -16}, {-190, -20}, {-310, -20}, {-310, -126}}, color = {95, 95, 95}));
    connect(bodyCylinder11.frame_a, revolute7.frame_b) annotation (
      Line(points = {{-144, -222}, {-144, -238.813}, {-174, -238.813}, {-174, -203}}, color = {95, 95, 95}));
    connect(bodyCylinder12.frame_a, revolute7.frame_b) annotation (
      Line(points = {{-146, -182}, {-146, -203.5}, {-174, -203.5}, {-174, -203}}, color = {95, 95, 95}));
    connect(fixedTranslation5.frame_b, revolute8.frame_a) annotation (
      Line(points = {{-300, -222}, {-300, -224}, {-286, -224}}, color = {95, 95, 95}));
    connect(fixedTranslation1.frame_b, fixedTranslation5.frame_a) annotation (
      Line(points = {{-188, 12}, {-320, 12}, {-320, -222}}, color = {95, 95, 95}));
    connect(bodyCylinder13.frame_a, revolute9.frame_b) annotation (
      Line(points = {{-132, 194}, {-132, 185.188}, {-144, 185.188}, {-144, 184.375}, {-164, 184.375}, {-164, 217}}, color = {95, 95, 95}));
    connect(bodyCylinder14.frame_a, revolute9.frame_b) annotation (
      Line(points = {{-128, 230}, {-128, 240.5}, {-164, 240.5}, {-164, 217}}, color = {95, 95, 95}));
    connect(revolute9.frame_a, fixedTranslation9.frame_b) annotation (
      Line(points={{-186,217},{-198,217},{-198,216}},        color = {95, 95, 95}));
    connect(revolute10.frame_b, fixedTranslation9.frame_a) annotation (
      Line(points = {{-254, 240}, {-229, 240}, {-229, 216}, {-218, 216}}));
    connect(fixedTranslation7.frame_b, revolute10.frame_a) annotation (
      Line(points = {{-276, 212}, {-274, 212}, {-274, 240}}, color = {95, 95, 95}));
    connect(revolute8.frame_b, fixedTranslation6.frame_a) annotation (
      Line(points = {{-266, -224}, {-240.5, -224}, {-240.5, -202}, {-236, -202}}));
    connect(revolute7.frame_a, fixedTranslation6.frame_b) annotation (
      Line(points = {{-196, -203}, {-216, -203}, {-216, -202}}, color = {95, 95, 95}));
    connect(revolute6.frame_a, fixedTranslation4.frame_b) annotation (
      Line(points = {{-192, -121}, {-208, -121}, {-208, -124}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_b, fixedTranslation10.frame_a) annotation (
      Line(points = {{-190, -16}, {-190, -32}, {-198, -32}, {-198, -48}}, color = {95, 95, 95}));
    connect(fixedTranslation8.frame_a, revolute5.frame_b) annotation (
      Line(points={{-298,-84},{-298,-114},{-248,-114},{-248,-125}}));
    connect(springDamperParallel.frame_b, fixedTranslation8.frame_b) annotation (
      Line(points = {{-264, -68}, {-252, -68}, {-252, -84}, {-278, -84}}, color = {95, 95, 95}));
    connect(springDamperParallel.frame_a, fixedTranslation10.frame_b) annotation (
      Line(points = {{-264, -48}, {-218, -48}}, color = {95, 95, 95}));
    connect(fixedTranslation12.frame_a, fixedTranslation1.frame_b) annotation (
      Line(points = {{-366, -10}, {-220, -10}, {-220, 12}, {-188, 12}}));
    connect(fixedTranslation12.frame_b, springDamperParallel1.frame_a) annotation (
      Line(points = {{-386, -10}, {-401, -10}, {-401, -72}}, color = {95, 95, 95}));
    connect(fixedTranslation11.frame_a, revolute8.frame_b) annotation (
      Line(points = {{-378, -182}, {-266, -182}, {-266, -224}}, color = {95, 95, 95}));
    connect(fixedTranslation11.frame_b, springDamperParallel1.frame_b) annotation (
      Line(points={{-398,-182},{-401,-182},{-401,-94}},        color = {95, 95, 95}));
    connect(fixedTranslation14.frame_a, fixedTranslation1.frame_b) annotation (
      Line(points = {{-340, 30}, {-188, 30}, {-188, 12}}, color = {95, 95, 95}));
    connect(fixedTranslation14.frame_b, springDamperParallel2.frame_a) annotation (
      Line(points = {{-360, 30}, {-410, 30}, {-410, 46}}, color = {95, 95, 95}));
    connect(revolute10.frame_b, fixedTranslation13.frame_a) annotation (
      Line(points = {{-254, 240}, {-224, 240}, {-224, 152}, {-250, 152}}));
    connect(fixedTranslation13.frame_b, springDamperParallel2.frame_b) annotation (
      Line(points = {{-250, 132}, {-410, 132}, {-410, 66}}));
    connect(bodyCylinder.frame_b, contact_sprocket.frame_b) annotation (
      Line(points={{-108,278},{16.2,278},{16.2,238.4}}));
    connect(revolute9.frame_b, contact.frame_b) annotation (
      Line(points={{-164,217},{-72,217},{-72,175.6},{-35.2,175.6}}));
    connect(contact.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{52,174.8},{204,174.8},{204,65.66}}, color = {95, 95, 95}, thickness = 0.5));
    connect(revolute4.frame_b, contact1.frame_b) annotation (
      Line(points={{-158,127},{-60,127},{-60,79.52},{-38.84,79.52}}, color = {95, 95, 95}));
    connect(contact1.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{44,78.76},{204,78.76},{204,65.66}},
                                                      thickness = 0.5));
    connect(revolute1.frame_b, contact2.frame_b) annotation (
      Line(points={{-156,55},{-58,55},{-58,-8.56},{-34.48,-8.56}}));
    connect(contact2.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{44,-9.28},{204,-9.28},{204,65.66}}, thickness = 0.5));
    connect(revolute2.frame_b, contact3.frame_b) annotation (
      Line(points={{-146,13},{-76,13},{-76,-102.24},{-61.92,-102.24}}, color = {95, 95, 95}));
    connect(contact3.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{34,-103.12},{204,-103.12},{204,65.66}},
                                                          color = {95, 95, 95}, thickness = 0.5));
    connect(revolute3.frame_b, contact4.frame_b) annotation (
      Line(points={{-162,-67},{-82,-67},{-82,-191.44},{-67.02,-191.44}}));
    connect(contact4.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{18,-192.22},{204,-192.22},{204,65.66}},
                                                          thickness = 0.5));
    connect(revolute6.frame_b, contact5.frame_b) annotation (
      Line(points={{-170,-121},{-106,-121},{-106,-286.24},{-69.92,-286.24}}, color = {95, 95, 95}));
    connect(contact5.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{26,-287.12},{204,-287.12},{204,65.66}},
                                                          color = {95, 95, 95}, thickness = 0.5));
    connect(revolute7.frame_b, contact6.frame_b) annotation (
      Line(points={{-174,-203},{-224,-203},{-224,-345.04},{-186.82,-345.04}}));
    connect(contact6.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{-80,-346.02},{204,-346.02},{204,65.66}},
                                                           thickness = 0.5));
    connect(contact_sprocket.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{38,238.2},{124,238.2},{124,65.66},{204,65.66}},thickness = 0.5));
    connect(fixedTranslation15.frame_a, fixedTranslation2.frame_b) annotation (
      Line(points = {{-562, 42}, {-586, 42}, {-586, -16}, {-190, -16}}));
    connect(fixedTranslation15.frame_b, revolute.frame_a) annotation (
      Line(points = {{-562, 62}, {-448, 62}, {-448, 270}, {-196, 270}}, color = {95, 95, 95}));
    connect(fixedTranslation.frame_a, fixedTranslation7.frame_a) annotation (
      Line(points = {{-208, 46}, {-320, 46}, {-320, 212}, {-296, 212}}, color = {95, 95, 95}));
    connect(fixedTranslation16.frame_b, fixedTranslation.frame_a) annotation (
      Line(points = {{-518, -44}, {-208, -44}, {-208, 46}}));
    connect(frame_a, fixedTranslation16.frame_a) annotation (
      Line(points = {{-664, -36}, {-538, -36}, {-538, -44}}));
    connect(cinghiafinal.frame_a, frame_b) annotation (
      Line(points={{338,61.64},{494,61.64},{494,64}},  thickness = 0.5));
    connect(damper.flange_b, revolute.axis) annotation (Line(points={{-178,316},
            {-162,316},{-162,280},{-186,280}}, color={0,0,0}));
    connect(damper.flange_a, revolute.support) annotation (Line(points={{-198,
            316},{-216,316},{-216,300},{-222,300},{-222,280},{-192,280}}, color=
           {0,0,0}));
    connect(prismatic.frame_b, revolute4.frame_a) annotation (Line(
        points={{-204,98},{-204,130},{-180,130},{-180,127}},
        color={95,95,95},
        thickness=0.5));
    connect(prismatic.frame_a, fixedTranslation.frame_a) annotation (Line(
        points={{-204,78},{-204,58},{-208,58},{-208,46}},
        color={95,95,95},
        thickness=0.5));
    connect(flange_a, revolute.axis) annotation (Line(points={{4,100},{-280,100},
            {-280,330},{-186,330},{-186,280}}, color={0,0,0}));
    connect(absoluteAngularVelocity1.w[3], filter.u) annotation (Line(points={{-107,
            334.333},{-88,334.333},{-88,334}}, color={0,0,127}));
    connect(absoluteAngularVelocity2.w[3], filter1.u) annotation (Line(points={{21,
            308.333},{40,308.333},{40,308}},     color={0,0,127}));
    connect(absoluteAngularVelocity1.frame_a, revolute.frame_b) annotation (
        Line(
        points={{-128,334},{-146,334},{-146,288},{-176,288},{-176,270}},
        color={95,95,95},
        thickness=0.5));
    connect(absoluteAngularVelocity2.frame_a, revolute9.frame_b) annotation (
        Line(
        points={{0,308},{-48,308},{-48,284},{-164,284},{-164,217}},
        color={95,95,95},
        thickness=0.5));
    connect(absoluteAngularVelocity3.w[3], filter2.u) annotation (Line(points={{57,
            268.333},{76,268.333},{76,268}},     color={0,0,127}));
    connect(absoluteAngularVelocity3.frame_a, revolute4.frame_b) annotation (
        Line(
        points={{36,268},{4,268},{4,250},{-158,250},{-158,127}},
        color={95,95,95},
        thickness=0.5));
    connect(absoluteAngularVelocity4.w[3], filter3.u) annotation (Line(points={{177,
            202.333},{196,202.333},{196,202}},      color={0,0,127}));
    connect(absoluteAngularVelocity4.frame_a, revolute1.frame_b) annotation (
        Line(
        points={{156,202},{76,202},{76,186},{-156,186},{-156,55}},
        color={95,95,95},
        thickness=0.5));
    connect(absoluteAngularVelocity5.w[3], filter4.u) annotation (Line(points={{337,
            -55.6667},{356,-55.6667},{356,-56}},      color={0,0,127}));
    connect(revolute2.frame_b, absoluteAngularVelocity5.frame_a) annotation (
        Line(
        points={{-146,13},{-108,13},{-108,8},{316,8},{316,-56}},
        color={95,95,95},
        thickness=0.5));
    connect(absoluteAngularVelocity6.w[3], filter5.u) annotation (Line(points={{275,
            -81.6667},{294,-81.6667},{294,-82}},      color={0,0,127}));
    connect(revolute3.frame_b, absoluteAngularVelocity6.frame_a) annotation (
        Line(
        points={{-162,-67},{-114,-67},{-114,-36},{254,-36},{254,-82}},
        color={95,95,95},
        thickness=0.5));
    connect(absoluteAngularVelocity7.w[3], filter6.u) annotation (Line(points={{275,
            -171.667},{294,-171.667},{294,-172}},      color={0,0,127}));
    connect(revolute6.frame_b, absoluteAngularVelocity7.frame_a) annotation (
        Line(
        points={{-170,-121},{-48,-121},{-48,-182},{254,-182},{254,-172}},
        color={95,95,95},
        thickness=0.5));
    connect(absoluteAngularVelocity8.w[3], filter7.u) annotation (Line(points={{287,
            -235.667},{306,-235.667},{306,-236}},      color={0,0,127}));
    connect(revolute7.frame_b, absoluteAngularVelocity8.frame_a) annotation (
        Line(
        points={{-174,-203},{-34,-203},{-34,-250},{266,-250},{266,-236}},
        color={95,95,95},
        thickness=0.5));
    connect(position.flange, prismatic.axis) annotation (Line(points={{-232,108},{
            -226,108},{-226,94},{-210,94},{-210,96}}, color={0,127,0}));
    connect(realExpression.y, position.s_ref) annotation (Line(points={{-257,82},{
            -256,82},{-256,108},{-254,108}}, color={0,0,127}));
    annotation (
      Icon(coordinateSystem(grid = {2, 0})),
      experiment(StartTime = 0, StopTime = 10, Tolerance = 0.0001, Interval = 0.02));
  end trackterrainfinal;

  model trackterrainfinal1
    parameter Real start[3] = {0, 0, 0};
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder(diameter = 0.16,
      length=0.105,
      lengthDirection={0,0,1},                                                                                               r = {0, 0, 0.0525}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-120, 276}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute(animation = false, useAxisFlange = true) annotation (
      Placement(visible = true, transformation(origin = {-186, 270}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder4(diameter = 0.088, length = 0.075, r = {0, 0, 0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-124, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute1(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-167, 55}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder1(diameter = 0.088, length = 0.075, r = {0, 0, -0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-122, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder2(diameter = 0.088, length = 0.075, r = {0, 0, 0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-125, -3}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute2(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-157, 11}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder5(diameter = 0.088, length = 0.075, r = {0, 0, -0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-129, -35}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder6(diameter = 0.088, length = 0.075, r = {0, 0, 0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-132, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute3(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-173, -67}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder7(diameter = 0.130, length = 0.075, r = {0, 0, 0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-124, 114}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder8(diameter = 0.130, length = 0.075, r = {0, 0, -0.0375},
      r_0(start={0.858,0.065,0.0525} + start),                                                                           useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-124, 160}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute4(animation = false) annotation (
      Placement(visible = true, transformation(origin={-167,127},    extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation(animation = true, r={-0.153,
          0.040,0})                                                                                                annotation (
      Placement(visible = true, transformation(origin={-200,46},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation1(r = {-0.350, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-188, 22}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation2(r = {-0.25, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-190, -6}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute5(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-261, -125}, extent = {{-13, 13}, {13, -13}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation3(animation = false, r = {0.250, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-300, -126}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation4(r = {-0.192, -0.192, 0}) annotation (
      Placement(visible = true, transformation(origin = {-218, -124}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder9(diameter = 0.105, length = 0.075, r = {0, 0, -0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-132, -110}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute6(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-181, -121}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder10(diameter = 0.105, length = 0.075, r = {0, 0, 0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-128, -148}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder3(diameter = 0.088, length = 0.075, r = {0, 0, -0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-122, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder11(diameter = 0.105, length = 0.075, r = {0, 0, 0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-134, -222}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute7(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-185, -203}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation5(animation = false, r = {0.281, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-310, -222}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder12(diameter = 0.105, length = 0.075, r = {0, 0, -0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-134, -182}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute8(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-276, -224}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation7(r = {-0.087, -0.029, 0}) annotation (
      Placement(visible = true, transformation(origin = {-288, 212}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder13(diameter = 0.105, length = 0.075, r = {0, 0, 0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-122, 194}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute9(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-175, 217}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder14(diameter = 0.105, length = 0.075, r = {0, 0, -0.0375}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-120, 230}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation9(r = {-0.071, -0.120, 0}) annotation (
      Placement(visible = true, transformation(origin = {-208, 216}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute10(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-246, 218}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation6(r = {-0.192, -0.192, 0}) annotation (
      Placement(visible = true, transformation(origin = {-226, -202}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation8(animation = false, r = {-0.134, -0.134, 0}) annotation (
      Placement(visible = true, transformation(origin = {-288, -84}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation10(animation = false, r = {0.064, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-208, -48}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel springDamperParallel(
      c=100000,
      d=1000,                                                                                             fixedRotationAtFrame_a = false, fixedRotationAtFrame_b = false,
      s_unstretched=0.15)                                                                                                                                                                       annotation (
      Placement(visible = true, transformation(origin = {-264, -58}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation11(animation = false, r = {-0.134, -0.134, 0}) annotation (
      Placement(visible = true, transformation(origin = {-388, -182}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation12(animation = false, r = {0.098, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-376, -10}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel springDamperParallel1(
      c=100000,
      d=1000,                                                                                              fixedRotationAtFrame_a = false, fixedRotationAtFrame_b = false,
      s_unstretched=0.14)                                                                                                                                                                        annotation (
      Placement(visible = true, transformation(origin = {-401, -83}, extent = {{-11, -11}, {11, 11}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel springDamperParallel2(
      c=100000,
      d=1000,                                                                                              fixedRotationAtFrame_a = false, fixedRotationAtFrame_b = false,
      s_unstretched=0.17)                                                                                                                                                                        annotation (
      Placement(visible = true, transformation(origin = {-410, 56}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation13(animation = false, r = {-0.044, -0.065, 0}) annotation (
      Placement(visible = true, transformation(origin = {-250, 142}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation14(animation = false, r = {0.278, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-350, 30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    AGILEX_BUNKER_PRO.contact_sprocket contact_sprocket annotation (Placement(
          visible=true, transformation(
          origin={28,238},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    AGILEX_BUNKER_PRO.cinghiaterrain cinghiafinal(start=start) annotation (
        Placement(visible=true, transformation(
          origin={271,67},
          extent={{-67,-67},{67,67}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact(R=0.0525) annotation (Placement(visible=
            true, transformation(
          origin={12,174},
          extent={{-40,-40},{40,40}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact1(R=0.065) annotation (Placement(visible=
            true, transformation(
          origin={6,78},
          extent={{-38,-38},{38,38}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact2(R=0.044) annotation (Placement(visible=
            true, transformation(
          origin={8,-10},
          extent={{-36,-36},{36,36}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact3(R=0.044) annotation (Placement(visible=
            true, transformation(
          origin={-10,-104},
          extent={{-44,-44},{44,44}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact4(R=0.044) annotation (Placement(visible=
            true, transformation(
          origin={-21,-193},
          extent={{-39,-39},{39,39}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact5(R=0.0525) annotation (Placement(visible=
            true, transformation(
          origin={-18,-288},
          extent={{-44,-44},{44,44}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact6(R=0.0525) annotation (Placement(visible=
            true, transformation(
          origin={-129,-347},
          extent={{-49,-49},{49,49}},
          rotation=0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation15(r={-0.105,
          -0.105,-0.0525})                                                                              annotation (
      Placement(visible = true, transformation(origin = {-562, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation16(animation = false, r={0.458,
          0.065,-0.0525})                                                                                                annotation (
      Placement(visible = true, transformation(origin = {-530, -46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a annotation (
      Placement(visible = true, transformation(origin = {-664, -36}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-100, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b[80] annotation (
      Placement(visible = true, transformation(origin = {494, 64}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {98, 0}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Components.Damper damper(d=7)
      annotation (Placement(transformation(extent={{-200,324},{-180,344}})));
    Modelica.Mechanics.MultiBody.Joints.Prismatic prismatic(useAxisFlange=true)
                       annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-202,96})));
    Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_a annotation (
        Placement(transformation(extent={{-6,86},{14,106}}), iconTransformation(
            extent={{-6,86},{14,106}})));
    Modelica.Mechanics.Translational.Sources.Position position(exact=true)
      annotation (Placement(transformation(extent={{-248,84},{-228,104}})));
    Modelica.Blocks.Sources.RealExpression realExpression
      annotation (Placement(transformation(extent={{-278,72},{-258,92}})));
  equation
    connect(revolute.frame_b, bodyCylinder.frame_a) annotation (
      Line(points = {{-176, 270}, {-140, 270}, {-140, 276}, {-130, 276}}));
    connect(bodyCylinder4.frame_a, revolute1.frame_b) annotation (
      Line(points = {{-134, 60}, {-134, 52.75}, {-138, 52.75}, {-138, 55.5}, {-156, 55.5}, {-156, 55}}, color = {95, 95, 95}));
    connect(bodyCylinder1.frame_a, revolute2.frame_b) annotation (
      Line(points = {{-132, 28}, {-146, 28}, {-146, 11}}, color = {95, 95, 95}));
    connect(bodyCylinder2.frame_a, revolute2.frame_b) annotation (
      Line(points = {{-134, -3}, {-134, -7.5}, {-146, -7.5}, {-146, 11}}, color = {95, 95, 95}));
    connect(bodyCylinder6.frame_a, revolute3.frame_b) annotation (
      Line(points = {{-142, -70}, {-142, -75.5}, {-162, -75.5}, {-162, -67}}, color = {95, 95, 95}));
    connect(bodyCylinder5.frame_a, revolute3.frame_b) annotation (
      Line(points = {{-140, -35}, {-140, -29.5}, {-166, -29.5}, {-166, -55.25}, {-162, -55.25}, {-162, -67}}, color = {95, 95, 95}));
    connect(bodyCylinder8.frame_a, revolute4.frame_b) annotation (
      Line(points={{-134,160},{-134,157.5},{-156,157.5},{-156,127}},          color = {95, 95, 95}));
    connect(bodyCylinder7.frame_a, revolute4.frame_b) annotation (
      Line(points={{-134,114},{-134,110.375},{-156,110.375},{-156,127}},          color = {95, 95, 95}));
    connect(fixedTranslation.frame_b, revolute1.frame_a) annotation (
      Line(points={{-190,46},{-190,55},{-178,55}}));
    connect(fixedTranslation.frame_b, fixedTranslation1.frame_a) annotation (
      Line(points={{-190,46},{-190,40},{-188,40},{-188,32}}));
    connect(fixedTranslation1.frame_b, revolute2.frame_a) annotation (
      Line(points = {{-188, 12}, {-188, 11}, {-168, 11}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_b, revolute3.frame_a) annotation (
      Line(points = {{-190, -16}, {-191.25, -16}, {-191.25, -32}, {-192.5, -32}, {-192.5, -40}, {-184, -40}, {-184, -67}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_a, fixedTranslation1.frame_b) annotation (
      Line(points = {{-190, 4}, {-189, 4}, {-189, 12}, {-188, 12}}));
    connect(fixedTranslation3.frame_b, revolute5.frame_a) annotation (
      Line(points = {{-290, -126}, {-303, -126}, {-303, -125}, {-274, -125}}, color = {95, 95, 95}));
    connect(revolute5.frame_b, fixedTranslation4.frame_a) annotation (
      Line(points = {{-248, -125}, {-248, -137.5}, {-228, -137.5}, {-228, -124}}));
    connect(bodyCylinder10.frame_a, revolute6.frame_b) annotation (
      Line(points = {{-138, -148}, {-138, -151.813}, {-146, -151.813}, {-146, -149.625}, {-170, -149.625}, {-170, -121}}, color = {95, 95, 95}));
    connect(bodyCylinder9.frame_a, revolute6.frame_b) annotation (
      Line(points = {{-142, -110}, {-142, -106.75}, {-160, -106.75}, {-160, -105.5}, {-170, -105.5}, {-170, -121}}, color = {95, 95, 95}));
    connect(bodyCylinder3.frame_a, revolute1.frame_b) annotation (
      Line(points = {{-132, 86}, {-156, 86}, {-156, 55}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_b, fixedTranslation3.frame_a) annotation (
      Line(points = {{-190, -16}, {-190, -20}, {-310, -20}, {-310, -126}}, color = {95, 95, 95}));
    connect(bodyCylinder11.frame_a, revolute7.frame_b) annotation (
      Line(points = {{-144, -222}, {-144, -238.813}, {-174, -238.813}, {-174, -203}}, color = {95, 95, 95}));
    connect(bodyCylinder12.frame_a, revolute7.frame_b) annotation (
      Line(points = {{-146, -182}, {-146, -203.5}, {-174, -203.5}, {-174, -203}}, color = {95, 95, 95}));
    connect(fixedTranslation5.frame_b, revolute8.frame_a) annotation (
      Line(points = {{-300, -222}, {-300, -224}, {-286, -224}}, color = {95, 95, 95}));
    connect(fixedTranslation1.frame_b, fixedTranslation5.frame_a) annotation (
      Line(points = {{-188, 12}, {-320, 12}, {-320, -222}}, color = {95, 95, 95}));
    connect(bodyCylinder13.frame_a, revolute9.frame_b) annotation (
      Line(points = {{-132, 194}, {-132, 185.188}, {-144, 185.188}, {-144, 184.375}, {-164, 184.375}, {-164, 217}}, color = {95, 95, 95}));
    connect(bodyCylinder14.frame_a, revolute9.frame_b) annotation (
      Line(points = {{-130, 230}, {-130, 240.5}, {-164, 240.5}, {-164, 217}}, color = {95, 95, 95}));
    connect(revolute9.frame_a, fixedTranslation9.frame_b) annotation (
      Line(points={{-186,217},{-198,217},{-198,216}},        color = {95, 95, 95}));
    connect(revolute10.frame_b, fixedTranslation9.frame_a) annotation (
      Line(points = {{-236, 218}, {-223, 218}, {-223, 216}, {-218, 216}}));
    connect(revolute8.frame_b, fixedTranslation6.frame_a) annotation (
      Line(points = {{-266, -224}, {-240.5, -224}, {-240.5, -202}, {-236, -202}}));
    connect(revolute7.frame_a, fixedTranslation6.frame_b) annotation (
      Line(points = {{-196, -203}, {-216, -203}, {-216, -202}}, color = {95, 95, 95}));
    connect(revolute6.frame_a, fixedTranslation4.frame_b) annotation (
      Line(points = {{-192, -121}, {-208, -121}, {-208, -124}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_b, fixedTranslation10.frame_a) annotation (
      Line(points = {{-190, -16}, {-190, -32}, {-198, -32}, {-198, -48}}, color = {95, 95, 95}));
    connect(fixedTranslation8.frame_a, revolute5.frame_b) annotation (
      Line(points={{-298,-84},{-298,-114},{-248,-114},{-248,-125}}));
    connect(springDamperParallel.frame_b, fixedTranslation8.frame_b) annotation (
      Line(points = {{-264, -68}, {-252, -68}, {-252, -84}, {-278, -84}}, color = {95, 95, 95}));
    connect(springDamperParallel.frame_a, fixedTranslation10.frame_b) annotation (
      Line(points = {{-264, -48}, {-218, -48}}, color = {95, 95, 95}));
    connect(fixedTranslation12.frame_a, fixedTranslation1.frame_b) annotation (
      Line(points = {{-366, -10}, {-220, -10}, {-220, 12}, {-188, 12}}));
    connect(fixedTranslation12.frame_b, springDamperParallel1.frame_a) annotation (
      Line(points = {{-386, -10}, {-401, -10}, {-401, -72}}, color = {95, 95, 95}));
    connect(fixedTranslation11.frame_a, revolute8.frame_b) annotation (
      Line(points = {{-378, -182}, {-266, -182}, {-266, -224}}, color = {95, 95, 95}));
    connect(fixedTranslation11.frame_b, springDamperParallel1.frame_b) annotation (
      Line(points = {{-398, -182}, {-401, -182}, {-401, -94}}, color = {95, 95, 95}));
    connect(fixedTranslation14.frame_a, fixedTranslation1.frame_b) annotation (
      Line(points={{-340,30},{-188,30},{-188,12}},        color = {95, 95, 95}));
    connect(fixedTranslation14.frame_b, springDamperParallel2.frame_a) annotation (
      Line(points = {{-360, 30}, {-410, 30}, {-410, 46}}, color = {95, 95, 95}));
    connect(revolute10.frame_b, fixedTranslation13.frame_a) annotation (
      Line(points = {{-236, 218}, {-224, 218}, {-224, 152}, {-250, 152}}));
    connect(fixedTranslation13.frame_b, springDamperParallel2.frame_b) annotation (
      Line(points = {{-250, 132}, {-410, 132}, {-410, 66}}));
    connect(bodyCylinder.frame_b, contact_sprocket.frame_b) annotation (
      Line(points={{-110,276},{16.2,276},{16.2,238.4}}));
    connect(revolute9.frame_b, contact.frame_b) annotation (
      Line(points={{-164,217},{-72,217},{-72,175.6},{-35.2,175.6}}));
    connect(contact.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{52,174.8},{204,174.8},{204,65.66}}, color = {95, 95, 95}, thickness = 0.5));
    connect(revolute4.frame_b, contact1.frame_b) annotation (
      Line(points={{-156,127},{-60,127},{-60,79.52},{-38.84,79.52}}, color = {95, 95, 95}));
    connect(contact1.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{44,78.76},{204,78.76},{204,65.66}},
                                                      thickness = 0.5));
    connect(revolute1.frame_b, contact2.frame_b) annotation (
      Line(points={{-156,55},{-58,55},{-58,-8.56},{-34.48,-8.56}}));
    connect(contact2.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{44,-9.28},{204,-9.28},{204,65.66}}, thickness = 0.5));
    connect(revolute2.frame_b, contact3.frame_b) annotation (
      Line(points={{-146,11},{-76,11},{-76,-102.24},{-61.92,-102.24}}, color = {95, 95, 95}));
    connect(contact3.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{34,-103.12},{204,-103.12},{204,65.66}},
                                                          color = {95, 95, 95}, thickness = 0.5));
    connect(revolute3.frame_b, contact4.frame_b) annotation (
      Line(points={{-162,-67},{-82,-67},{-82,-191.44},{-67.02,-191.44}}));
    connect(contact4.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{18,-192.22},{204,-192.22},{204,65.66}},
                                                          thickness = 0.5));
    connect(revolute6.frame_b, contact5.frame_b) annotation (
      Line(points={{-170,-121},{-106,-121},{-106,-286.24},{-69.92,-286.24}}, color = {95, 95, 95}));
    connect(contact5.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{26,-287.12},{204,-287.12},{204,65.66}},
                                                          color = {95, 95, 95}, thickness = 0.5));
    connect(revolute7.frame_b, contact6.frame_b) annotation (
      Line(points={{-174,-203},{-224,-203},{-224,-345.04},{-186.82,-345.04}}));
    connect(contact6.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{-80,-346.02},{204,-346.02},{204,65.66}},
                                                           thickness = 0.5));
    connect(contact_sprocket.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{38,238.2},{124,238.2},{124,65.66},{204,65.66}},thickness = 0.5));
  connect(fixedTranslation15.frame_a, fixedTranslation2.frame_b) annotation (
      Line(points = {{-562, 42}, {-586, 42}, {-586, -16}, {-190, -16}}));
  connect(fixedTranslation15.frame_b, revolute.frame_a) annotation (
      Line(points = {{-562, 62}, {-448, 62}, {-448, 270}, {-196, 270}}, color = {95, 95, 95}));
    connect(fixedTranslation16.frame_b, fixedTranslation.frame_a) annotation (
      Line(points={{-520,-46},{-210,-46},{-210,46}}));
    connect(frame_a, fixedTranslation16.frame_a) annotation (
      Line(points = {{-664, -36}, {-540, -36}, {-540, -46}}));
    connect(cinghiafinal.frame_a, frame_b) annotation (
      Line(points={{338,61.64},{494,61.64},{494,64}},  thickness = 0.5));
    connect(fixedTranslation7.frame_b, revolute10.frame_a) annotation (
      Line(points = {{-278, 212}, {-256, 212}, {-256, 218}}, color = {95, 95, 95}));
    connect(damper.flange_b, revolute.axis) annotation (Line(points={{-180,334},
            {-180,307},{-186,307},{-186,280}}, color={0,0,0}));
    connect(damper.flange_a, revolute.support) annotation (Line(points={{-200,
            334},{-200,307},{-192,307},{-192,280}}, color={0,0,0}));
    connect(fixedTranslation7.frame_a, fixedTranslation.frame_a) annotation (
        Line(
        points={{-298,212},{-358,212},{-358,46},{-210,46}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute4.frame_a, prismatic.frame_b) annotation (Line(
        points={{-178,127},{-202,127},{-202,106}},
        color={95,95,95},
        thickness=0.5));
    connect(prismatic.frame_a, fixedTranslation.frame_a) annotation (Line(
        points={{-202,86},{-202,46},{-210,46}},
        color={95,95,95},
        thickness=0.5));
    connect(flange_a, revolute.axis) annotation (Line(points={{4,96},{-424,96},
            {-424,350},{-186,350},{-186,280}}, color={0,0,0}));
    connect(position.flange, prismatic.axis) annotation (Line(points={{-228,94},{-230,
            94},{-230,104},{-208,104}}, color={0,127,0}));
    connect(realExpression.y, position.s_ref) annotation (Line(points={{-257,82},{
            -252,82},{-252,90},{-250,90},{-250,94}}, color={0,0,127}));
    annotation (
      Icon(coordinateSystem(grid = {2, 0})),
      experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-4, Interval = 0.02));
  end trackterrainfinal1;

  model completevehicle
  extends Modelica.Icons.Example;
    AGILEX_BUNKER_PRO.trackterrainfinal trackterrainfinal(start={2.7,0.15,0})
      annotation (Placement(visible=true, transformation(
          origin={60,18},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact_t contact_t annotation (Placement(visible=true,
          transformation(
          origin={110,18},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    AGILEX_BUNKER_PRO.trackterrainfinal1 trackterrainfinal1(start={2.7,0.15,-0.90})
      annotation (Placement(visible=true, transformation(
          origin={-50,18},
          extent={{10,-10},{-10,10}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact_t contact_t1 annotation (Placement(visible=true,
          transformation(
          origin={-98,18},
          extent={{10,-10},{-10,10}},
          rotation=0)));
    inner AGILEX_BUNKER_PRO.Surface Surface(
      nu=4,
      nv=4,
      multiColoredSurface=false,
      color={0,255,0},
      FlatSurface=true,
      nu1=4,
      nv1=4,
      LengthX=8,
      LengthZ=8,
      d=sum(contact_t.contact_terrain[i].dd for i in 1:25) + sum(contact_t1.contact_terrain[
          i].dd for i in 1:25)) annotation (Placement(visible=true,
          transformation(
          origin={106,69},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    inner Modelica.Mechanics.MultiBody.World world annotation (
      Placement(visible = true, transformation(origin={-64,-51},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation(animation = false, r = {0, 0, 0.80}) annotation (
      Placement(visible = true, transformation(origin = {-2, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Sensors.AbsoluteAngles absoluteAngles annotation (
      Placement(visible = true, transformation(origin = {20, -27}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp(duration = 3,
      height=100,                                               startTime = 1) annotation (
      Placement(visible = true, transformation(origin = {-76, 63}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp1(duration = 3,
      height=100,                                                startTime = 10) annotation (
      Placement(visible = true, transformation(origin = {-6, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.Add add(k1 = -1) annotation (
      Placement(visible = true, transformation(origin={30,57},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Sources.Torque torque
      annotation (Placement(transformation(extent={{-60,36},{-40,57}})));
    Modelica.Mechanics.Rotational.Sources.Torque torque1
      annotation (Placement(transformation(extent={{56,45},{76,66}})));
  equation
    connect(trackterrainfinal.frame_b, contact_t.frame_b) annotation (
      Line(points={{69.8,18},{84,18},{84,18.4},{98.2,18.4}},
                                          color = {95, 95, 95}, thickness = 0.5));
    connect(contact_t1.frame_b, trackterrainfinal1.frame_b) annotation (
      Line(points={{-86.2,18.4},{-73.8,18.4},{-73.8,18},{-59.8,18}},
                                            color = {95, 95, 95}, thickness = 0.5));
    connect(trackterrainfinal1.frame_a, fixedTranslation.frame_a) annotation (
      Line(points={{-40,17.8},{-26,17.8},{-26,18},{-12,18}},
                                            color = {95, 95, 95}));
    connect(fixedTranslation.frame_b, trackterrainfinal.frame_a) annotation (
      Line(points={{8,18},{30,18},{30,17.8},{50,17.8}},
                                         color = {95, 95, 95}));
    connect(fixedTranslation.frame_b, absoluteAngles.frame_a) annotation (
      Line(points = {{8, 18}, {8, -4.5}, {10, -4.5}, {10, -27}}));
    connect(ramp1.y, add.u1) annotation (
      Line(points={{5,66},{10,66},{10,63},{18,63}},color = {0, 0, 127}));
    connect(ramp.y, add.u2) annotation (
      Line(points={{-65,63},{-22,63},{-22,51},{18,51}},          color = {0, 0, 127}));
    connect(ramp.y, torque.tau) annotation (Line(points={{-65,63},{-62,63},{-62,
            46.5}}, color={0,0,127}));
    connect(torque.flange, trackterrainfinal1.flange_a) annotation (Line(points=
           {{-40,46.5},{-40,27.6},{-50.4,27.6}}, color={0,0,0}));
    connect(add.y, torque1.tau) annotation (Line(points={{41,57},{48,57},{48,
            55.5},{54,55.5}}, color={0,0,127}));
    connect(torque1.flange, trackterrainfinal.flange_a) annotation (Line(points=
           {{76,55.5},{84,55.5},{84,51},{92,51},{92,28},{60.4,28}}, color={0,0,
            0}));
    annotation (
      Icon(coordinateSystem(grid = {2, 0})),
      Diagram(coordinateSystem(grid = {2, 3}), graphics={Polygon(points={{108,
                63},{108,63}}, lineColor={28,108,200})}),
      experiment(
        StopTime=20,
        Interval=0.008,
        Tolerance=0.001,
        __Dymola_Algorithm="Dassl"));
  end completevehicle;

  model completevehiclevel
  extends Modelica.Icons.Example;
    AGILEX_BUNKER_PRO.trackterrainfinal trackterrainfinal(start={2.7,0.15,0})
      annotation (Placement(visible=true, transformation(
          origin={60,18},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact_t contact_t annotation (Placement(visible=true,
          transformation(
          origin={110,18},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    AGILEX_BUNKER_PRO.trackterrainfinal1 trackterrainfinal1(start={2.7,0.15,-0.90})
      annotation (Placement(visible=true, transformation(
          origin={-50,18},
          extent={{10,-10},{-10,10}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact_t contact_t1 annotation (Placement(visible=true,
          transformation(
          origin={-98,18},
          extent={{10,-10},{-10,10}},
          rotation=0)));
    inner AGILEX_BUNKER_PRO.Surface Surface(
      nu=4,
      nv=4,
      multiColoredSurface=false,
      color={0,255,0},
      FlatSurface=true,
      nu1=4,
      nv1=4,
      LengthX=8,
      LengthZ=8,
      d=sum(contact_t.contact_terrain[i].dd for i in 1:25) + sum(contact_t1.contact_terrain[
          i].dd for i in 1:25)) annotation (Placement(visible=true,
          transformation(
          origin={106,69},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    inner Modelica.Mechanics.MultiBody.World world annotation (
      Placement(visible = true, transformation(origin={-64,-51},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation(animation = false, r = {0, 0, 0.80}) annotation (
      Placement(visible = true, transformation(origin = {-2, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Sensors.AbsoluteAngles absoluteAngles annotation (
      Placement(visible = true, transformation(origin = {20, -27}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp(duration = 3,
      height=7,                                                 startTime = 1) annotation (
      Placement(visible = true, transformation(origin={-118,75},   extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp1(duration = 3,
      height=7,
      startTime=7)                                                               annotation (
      Placement(visible = true, transformation(origin = {-6, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.Add add(k1 = -1) annotation (
      Placement(visible = true, transformation(origin={36,60},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Sources.Speed speed
      annotation (Placement(transformation(extent={{60,50},{80,70}})));
    Modelica.Mechanics.Rotational.Sources.Speed speed1
      annotation (Placement(transformation(extent={{-78,43},{-58,63}})));
  equation
    connect(trackterrainfinal.frame_b, contact_t.frame_b) annotation (
      Line(points={{69.8,18},{84,18},{84,18.4},{98.2,18.4}},
                                          color = {95, 95, 95}, thickness = 0.5));
    connect(contact_t1.frame_b, trackterrainfinal1.frame_b) annotation (
      Line(points={{-86.2,18.4},{-73.8,18.4},{-73.8,18},{-59.8,18}},
                                            color = {95, 95, 95}, thickness = 0.5));
    connect(trackterrainfinal1.frame_a, fixedTranslation.frame_a) annotation (
      Line(points={{-40,17.8},{-26,17.8},{-26,18},{-12,18}},
                                            color = {95, 95, 95}));
    connect(fixedTranslation.frame_b, trackterrainfinal.frame_a) annotation (
      Line(points={{8,18},{30,18},{30,17.8},{50,17.8}},
                                         color = {95, 95, 95}));
    connect(fixedTranslation.frame_b, absoluteAngles.frame_a) annotation (
      Line(points = {{8, 18}, {8, -4.5}, {10, -4.5}, {10, -27}}));
    connect(ramp1.y, add.u1) annotation (
      Line(points={{5,66},{24,66}},                color = {0, 0, 127}));
    connect(ramp.y, add.u2) annotation (
      Line(points={{-107,75},{-22,75},{-22,54},{24,54}},         color = {0, 0, 127}));
    connect(ramp.y, speed1.w_ref) annotation (Line(points={{-107,75},{-96,75},{
            -96,60},{-80,60},{-80,53}}, color={0,0,127}));
    connect(speed1.flange, trackterrainfinal1.flange_a) annotation (Line(points=
           {{-58,53},{-48,53},{-48,27.6},{-50.4,27.6}}, color={0,0,0}));
    connect(add.y, speed.w_ref)
      annotation (Line(points={{47,60},{58,60}}, color={0,0,127}));
    connect(speed.flange, trackterrainfinal.flange_a) annotation (Line(points={
            {80,60},{94,60},{94,33},{60.4,33},{60.4,28}}, color={0,0,0}));
    annotation (
      Icon(coordinateSystem(grid = {2, 0})),
      Diagram(coordinateSystem(grid = {2, 3}), graphics={Polygon(points={{108,
                63},{108,63}}, lineColor={28,108,200})}),
      experiment(
        StopTime=20,
        Interval=0.008,
        Tolerance=0.001,
        __Dymola_Algorithm="Dassl"));
  end completevehiclevel;

  model completevehiclewithchasis
  extends Modelica.Icons.Example;
    AGILEX_BUNKER_PRO.trackterrainfinal trackterrainfinal(start={2.7,0.15,0})
      annotation (Placement(visible=true, transformation(
          origin={62,18},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact_t contact_t annotation (Placement(visible=true,
          transformation(
          origin={110,18},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    AGILEX_BUNKER_PRO.trackterrainfinal1 trackterrainfinal1(start={2.7,0.15,-0.905})
      annotation (Placement(visible=true, transformation(
          origin={-50,18},
          extent={{10,-10},{-10,10}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact_t contact_t1 annotation (Placement(visible=true,
          transformation(
          origin={-98,18},
          extent={{10,-10},{-10,10}},
          rotation=0)));
    inner AGILEX_BUNKER_PRO.Surface Surface(
      nu=2,
      nv=2,
      multiColoredSurface=false,
      color={0,255,0},
      FlatSurface=true,
      nu1=2,
      nv1=2,
      LengthX=8,
      LengthZ=8,
      d=sum(contact_t.contact_terrain[i].dd for i in 1:25) + sum(contact_t1.contact_terrain[
          i].dd for i in 1:25)) annotation (Placement(visible=true,
          transformation(
          origin={106,69},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    inner Modelica.Mechanics.MultiBody.World world annotation (
      Placement(visible = true, transformation(origin={-62,-45},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Sensors.AbsoluteAngles absoluteAngles annotation (
      Placement(visible = true, transformation(origin = {20, -27}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp(
      duration=0,
      height=250,                                               startTime = 1) annotation (
      Placement(visible = true, transformation(origin = {-76, 63}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp1(
      duration=0,
      height=500,
      startTime=1)                                                               annotation (
      Placement(visible = true, transformation(origin = {-6, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.Add add(k1 = -1) annotation (
      Placement(visible = true, transformation(origin={30,57},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Sources.Torque torque
      annotation (Placement(transformation(extent={{-60,36},{-40,57}})));
    Modelica.Mechanics.Rotational.Sources.Torque torque1
      annotation (Placement(transformation(extent={{56,45},{76,66}})));
    Modelica.Mechanics.MultiBody.Parts.BodyBox bodyBox(
      r={0,0,0.8},
      length=0.8,
      width=0.1,
      height=1,
      density=100)
      annotation (Placement(transformation(extent={{-12,8},{8,28}})));
  equation
    connect(trackterrainfinal.frame_b, contact_t.frame_b) annotation (
      Line(points={{71.8,18},{84,18},{84,18.4},{98.2,18.4}},
                                          color = {95, 95, 95}, thickness = 0.5));
    connect(contact_t1.frame_b, trackterrainfinal1.frame_b) annotation (
      Line(points={{-86.2,18.4},{-73.8,18.4},{-73.8,18},{-59.8,18}},
                                            color = {95, 95, 95}, thickness = 0.5));
    connect(ramp1.y, add.u1) annotation (
      Line(points={{5,66},{10,66},{10,63},{18,63}},color = {0, 0, 127}));
    connect(ramp.y, add.u2) annotation (
      Line(points={{-65,63},{-22,63},{-22,51},{18,51}},          color = {0, 0, 127}));
    connect(ramp.y, torque.tau) annotation (Line(points={{-65,63},{-62,63},{-62,
            46.5}}, color={0,0,127}));
    connect(torque.flange, trackterrainfinal1.flange_a) annotation (Line(points=
           {{-40,46.5},{-40,27.6},{-50.4,27.6}}, color={0,0,0}));
    connect(add.y, torque1.tau) annotation (Line(points={{41,57},{48,57},{48,
            55.5},{54,55.5}}, color={0,0,127}));
    connect(torque1.flange, trackterrainfinal.flange_a) annotation (Line(points={{76,55.5},
            {84,55.5},{84,51},{92,51},{92,28},{62.4,28}},           color={0,0,
            0}));
    connect(trackterrainfinal1.frame_a, bodyBox.frame_a) annotation (Line(
        points={{-40,17.8},{-32,17.8},{-32,16},{-12,16},{-12,18}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyBox.frame_b, trackterrainfinal.frame_a) annotation (Line(
        points={{8,18},{30,18},{30,17.8},{52,17.8}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyBox.frame_b, absoluteAngles.frame_a) annotation (Line(
        points={{8,18},{8,-27},{10,-27}},
        color={95,95,95},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=20,
        Interval=0.008,
        Tolerance=0.001,
        __Dymola_Algorithm="Dassl"));
  end completevehiclewithchasis;

  model completevehiclewithchasisvel
  extends Modelica.Icons.Example;
    AGILEX_BUNKER_PRO.trackterrainfinal trackterrainfinal(start={2.7,0.15,0})
      annotation (Placement(visible=true, transformation(
          origin={60,18},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact_t contact_t annotation (Placement(visible=true,
          transformation(
          origin={110,18},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    AGILEX_BUNKER_PRO.trackterrainfinal1 trackterrainfinal1(start={2.7,0.15,-0.905})
      annotation (Placement(visible=true, transformation(
          origin={-50,18},
          extent={{10,-10},{-10,10}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact_t contact_t1 annotation (Placement(visible=true,
          transformation(
          origin={-98,18},
          extent={{10,-10},{-10,10}},
          rotation=0)));
    inner AGILEX_BUNKER_PRO.Surface Surface(
      nu=2,
      nv=2,
      multiColoredSurface=false,
      color={0,255,0},
      FlatSurface=true,
      nu1=2,
      nv1=2,
      LengthX=8,
      LengthZ=8,
      d=sum(contact_t.contact_terrain[i].dd for i in 1:25) + sum(contact_t1.contact_terrain[
          i].dd for i in 1:25)) annotation (Placement(visible=true,
          transformation(
          origin={106,69},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    inner Modelica.Mechanics.MultiBody.World world annotation (
      Placement(visible = true, transformation(origin={-76,-45},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Sensors.AbsoluteAngles absoluteAngles annotation (
      Placement(visible = true, transformation(origin = {20, -27}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp(duration = 3,
      height=10,                                                startTime = 1) annotation (
      Placement(visible = true, transformation(origin={-118,75},   extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp1(duration = 3,
      height=20,
      startTime=1)                                                               annotation (
      Placement(visible = true, transformation(origin = {-6, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Math.Add add(k1 = -1) annotation (
      Placement(visible = true, transformation(origin={36,60},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Sources.Speed speed
      annotation (Placement(transformation(extent={{60,50},{80,70}})));
    Modelica.Mechanics.Rotational.Sources.Speed speed1
      annotation (Placement(transformation(extent={{-78,43},{-58,63}})));
    Modelica.Mechanics.MultiBody.Parts.BodyBox bodyBox(
      r={0,0,0.8},
      length=0.8,
      width=0.1,
      height=1,
      density=100)
      annotation (Placement(transformation(extent={{-6,8},{14,28}})));
  equation
    connect(trackterrainfinal.frame_b, contact_t.frame_b) annotation (
      Line(points={{69.8,18},{84,18},{84,18.4},{98.2,18.4}},
                                          color = {95, 95, 95}, thickness = 0.5));
    connect(contact_t1.frame_b, trackterrainfinal1.frame_b) annotation (
      Line(points={{-86.2,18.4},{-73.8,18.4},{-73.8,18},{-59.8,18}},
                                            color = {95, 95, 95}, thickness = 0.5));
    connect(ramp1.y, add.u1) annotation (
      Line(points={{5,66},{24,66}},                color = {0, 0, 127}));
    connect(ramp.y, add.u2) annotation (
      Line(points={{-107,75},{-22,75},{-22,54},{24,54}},         color = {0, 0, 127}));
    connect(ramp.y, speed1.w_ref) annotation (Line(points={{-107,75},{-96,75},{
            -96,60},{-80,60},{-80,53}}, color={0,0,127}));
    connect(speed1.flange, trackterrainfinal1.flange_a) annotation (Line(points=
           {{-58,53},{-48,53},{-48,27.6},{-50.4,27.6}}, color={0,0,0}));
    connect(add.y, speed.w_ref)
      annotation (Line(points={{47,60},{58,60}}, color={0,0,127}));
    connect(speed.flange, trackterrainfinal.flange_a) annotation (Line(points={
            {80,60},{94,60},{94,33},{60.4,33},{60.4,28}}, color={0,0,0}));
    connect(bodyBox.frame_b, trackterrainfinal.frame_a) annotation (Line(
        points={{14,18},{26,18},{26,14},{50,14},{50,17.8}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyBox.frame_a, trackterrainfinal1.frame_a) annotation (Line(
        points={{-6,18},{-23,18},{-23,17.8},{-40,17.8}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyBox.frame_b, absoluteAngles.frame_a) annotation (Line(
        points={{14,18},{4,18},{4,-27},{10,-27}},
        color={95,95,95},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=20,
        Interval=0.008,
        Tolerance=0.001,
        __Dymola_Algorithm="Dassl"));
  end completevehiclewithchasisvel;

  model Surface_base
     extends Modelica.Mechanics.MultiBody.Visualizers.Advanced.Surface;
    //"implemetation of the surface class based on interoplation between four elevation points with eight paritial derivatives."
    //import SI = Modelica.Units.SI;
    //import MB = Modelica.Mechanics.MultiBody;
    parameter Boolean FlatSurface = true;
    //"If true simpler equations for the functions can be used, false enables uneven surfaces";
    parameter Boolean DefaultCharSurface = true;
    parameter Boolean visSurface = true;
    //"= true if the surfrace shall be shown in the animation window";
    parameter Integer nu1(final min = 2) = 4;
    //"Number of grid points in x direction";
    parameter Integer nv1(final min = 2) = 4;
    //"Number of grid points in z direction";
    parameter Integer IP(final min = 0) = 5;
    //"Number of interpolation points between the grid points (just for visualization)";
    parameter Real PNG = 1;
    // "Filename of the that shall be used as texture";
    //parameter MB.Types.Color Color = {192, 192, 192};
    // "Surface color (mixed with texture)";
    //parameter MB.Types.SpecularCoefficient SpecularCoefficient = 0.1;
    // "Specular coefficient of the road surface without texture";
    parameter Modelica.Units.SI.Length LengthX = 50;
    //"Length of the surface area in x direction";
    parameter Modelica.Units.SI.Length LengthZ = 50;
    //"Length of the surface area in z direction";
    parameter Modelica.Units.SI.Length OffsetX = -LengthX / 2;
    // "Offset in x direction";
    parameter Modelica.Units.SI.Length OffsetZ = -LengthZ / 2;
    // "Offset in z direction";
    parameter Modelica.Units.SI.Position yel[nu1, nv1] = zeros(nu1, nv1);
    //"Matrix with the y (elevation) values at the grid points"
    parameter Real dy_dx[nu1, nv1] = zeros(nu1, nv1);
    //"Matrix with the slope of the grid points in x direction"
    parameter Real dy_dz[nu1, nv1] = zeros(nu1, nv1);
    //"Matrix with the slope of the grid points in y direction"
    parameter Real mu[nu1 - 1, nv1 - 1] = ones(nu1 - 1, nv1 - 1)*0.7;
    //"Matrix with the friction coefficient for every rectangle";
    parameter Real kc[nu1 - 1, nv1 - 1] = zeros(nu1 - 1, nv1 - 1);
    parameter Real kf[nu1 - 1, nv1 - 1] = ones(nu1 - 1, nv1 - 1)*1e10;
    parameter Real kv[nu1 - 1, nv1 - 1] = ones(nu1 - 1, nv1 - 1)*1e5;
    parameter Real n[nu1 - 1, nv1 - 1] = ones(nu1 - 1, nv1 - 1);
    parameter Modelica.Units.SI.ShearStress c[nu1 - 1, nv1 - 1] = ones(nu1 - 1, nv1 - 1)*2.2539e6;
    parameter Real tanf[nu1 - 1, nv1 - 1] = ones(nu1 - 1, nv1 - 1)*Modelica.Math.tan((44.5/180)*Modelica.Constants.pi);
    parameter Modelica.Units.SI.Length K[nu1 - 1, nv1 - 1] = ones(nu1 - 1, nv1 - 1)*0.0002;
    parameter Real gammas[nu1 - 1, nv1 - 1] = ones(nu1 - 1, nv1 - 1)*14126.4;
    final parameter Integer nuIP = nu1 + (nu1 - 1) * IP;
    // "Number of the grid values in x direction with interpolation";
    final parameter Integer nvIP = nv1 + (nv1 - 1) * IP;
    // "Number of the grid values in y direction with interpolation";
    function get_eN = get_eN_protected(FlatSurface = FlatSurface, yMatrix = yel, dy_dxMatrix = dy_dx, dy_dzMatrix = dy_dz, nu = nu1, nv = nv1, LengthX = LengthX, LengthZ = LengthZ, OffsetX = OffsetX, OffsetZ = OffsetZ) annotation (
      Documentation(info = "<html>
  
  See the documentation <a href=\"Modelica://WheelsAndTires.Environment.Surface\">Surface</a> for more information.
  
  </html>"));
    replaceable function get_elevation = get_elevation_protected(FlatSurface = FlatSurface, yMatrix = yel, dy_dxMatrix = dy_dx, dy_dzMatrix = dy_dz, nu = nu1, nv = nv1, LengthX = LengthX, LengthZ = LengthZ, OffsetX = OffsetX, OffsetZ = OffsetZ) constrainedby
      get_elevation_base                                                                                                                                                                                                         annotation (
      Documentation(info = "<html>
  
  See the documentation <a href=\"Modelica://WheelsAndTires.Environment.Surface\">Surface</a> for more information.
  
  </html>"));
    replaceable function get_surface_char = get_surface_char_protected(DefaultCharSurface = DefaultCharSurface, muMatrix = mu, kcMatrix = kc, kfMatrix = kf, kvMatrix = kv, nMatrix = n, cMatrix = c, tanfMatrix = tanf, KMatrix = K, gammasMatrix = gammas, nu = nu1, nv = nv1, LengthX = LengthX, LengthZ = LengthZ, OffsetX = OffsetX, OffsetZ = OffsetZ) constrainedby get_surface_char_base annotation (
      Documentation(info = "<html>
  
  See the documentation <a href=\"Modelica://WheelsAndTires.Environment.Surface\">Surface</a> for more information.
  
  </html>"));
    redeclare replaceable function surfaceCharacteristic = SurfaceFunction(x_min=OffsetX + 1e-6, x_max=LengthX+OffsetX - 1e-6, z_min=OffsetZ + 1e-6, z_max=LengthZ+OffsetZ - 1e-6) constrainedby
      Modelica.Mechanics.MultiBody.Interfaces.partialSurfaceCharacteristic;

    function get_rectangle
      //"finds the rectangle which calculates the elevation and normal vector for the current position of the tire."
      input Modelica.Units.SI.Position x;
      //"acutal position to find active recthangle and relative position in x direction";
      input Modelica.Units.SI.Position z;
      //"acutal position to find active recthangle and relative position in z direction";
      input Integer nu;
      // "Number of grid points in x direction";
      input Integer nv;
      // "Number of grid points in z direction";
      input Modelica.Units.SI.Length LengthX;
      // "Length of the surface area in x direction";
      input Modelica.Units.SI.Length LengthZ;
      // "Length of the surface area in z direction";
      input Modelica.Units.SI.Distance OffsetX;
      // "Offset of the surface area in x direction";
      input Modelica.Units.SI.Distance OffsetZ;
      // "Offset of the surface area in z direction";
      output Integer rectX;
      // "index of the active rectange in x direction";
      output Integer rectZ;
      // "index of the active rectange in z direction";
      output Boolean inArea;
      //"true: tire is in the defined surface area, false: it is not.";
      output Modelica.Units.SI.Distance localX;
      //"relative position in the active rectangle in x direction.";
      output Modelica.Units.SI.Distance localZ;
      //"relative position in the active rectangle in z direction.";
    protected
      Modelica.Units.SI.Distance rectangleX;
      Modelica.Units.SI.Distance rectangleZ;
    algorithm
      rectX := integer(ceil((x - OffsetX) / LengthX * (nu - 1)));
      rectZ := integer(ceil((z - OffsetZ) / LengthZ * (nv - 1)));
      if rectX > nu - 1 or rectX < 1 or rectZ > nv - 1 or rectZ < 1 then
        inArea := false;
      else
        inArea := true;
      end if;
      rectangleX := LengthX / (nu - 1);
      rectangleZ := LengthZ / (nv - 1);
      // getting local position (in the rectangle the tire is in at the time of calling)
      localX := mod(x - LengthX - OffsetX, rectangleX);
      localZ := mod(z - LengthZ - OffsetZ, rectangleZ);
    end get_rectangle;

  protected
  partial function get_elevation_base
      input Modelica.Units.SI.Position x;
      //"actual position to get the elevation of the surface (y-coordinate) in x direction"
      input Modelica.Units.SI.Position z;
      //"actual position to get the elevation of the surface (y-coordinate) in z direction"
      output Modelica.Units.SI.Position elevation;
  end get_elevation_base;

    partial function get_surface_char_base
      input Modelica.Units.SI.Position x;
      //"actual position to get the friction coefficient of the surface in x direction"
      input Modelica.Units.SI.Position z;
      //"actual position to get the friction coefficient of the surface in z direction"
      input Boolean DefaultCharSurface;
      input Real muMatrix[:, :];
      input Real kcMatrix[:, :];
      input Real kfMatrix[:, :];
      input Real kvMatrix[:, :];
      input Real nMatrix[:, :];
      input Modelica.Units.SI.ShearStress cMatrix[:, :];
      input Real tanfMatrix[:, :];
      input Modelica.Units.SI.Length KMatrix[:, :];
      input Real gammasMatrix[:, :];
      output Real mu;
      //"the friction coefficient";
      output Real kc;
      output Real kf;
      output Real kv;
      output Real n;
      output Modelica.Units.SI.ShearStress c;
      output Real tanf;
      output Modelica.Units.SI.Length K;
      output Real gammas;
    end get_surface_char_base;

    function get_eN_protected
      input Modelica.Units.SI.Position x;
      //"actual position to get the normal vector eN in x direction."
      input Modelica.Units.SI.Position z;
      //"actual position to get the normal vector eN in z direction."
      input Boolean FlatSurface;
      input Modelica.Units.SI.Distance yMatrix[:, :];
      input Real dy_dxMatrix[:, :];
      input Real dy_dzMatrix[:, :];
      input Integer nu;
      input Integer nv;
      input Modelica.Units.SI.Length LengthX;
      input Modelica.Units.SI.Length LengthZ;
      input Modelica.Units.SI.Distance OffsetX;
      input Modelica.Units.SI.Distance OffsetZ;
      output Real eN[3];
      // "normal vector on the surface.";
    protected
      Integer xIndex;
      Integer zIndex;
      Modelica.Units.SI.Distance localX;
      Modelica.Units.SI.Distance localZ;
      Boolean inArea;
      Modelica.Units.SI.Distance rectangleY[2, 2];
      Real rectangleDy_dx[2, 2];
      Real rectangleDy_dz[2, 2];
      Real X[4];
      Real Z[4];
      Real G[4, 4];
      constant Real H[4, 4] = [2, -3, 0, 1; -2, 3, 0, 0; 1, -2, 1, 0; 1, -1, 0, 0];
      constant Real dH[4, 4] = [0, 6, -6, 0; 0, -6, 6, 0; 0, 3, -4, 1; 0, 3, -2, 0];
      Real dy_dx;
      Real dy_dz;
      Modelica.Units.SI.Distance factorX;
      Modelica.Units.SI.Distance factorZ;
      Modelica.Units.SI.Length localXNorm;
      Modelica.Units.SI.Length localZNorm;
      Modelica.Units.SI.Length rectangleX;
      Modelica.Units.SI.Length rectangleZ;
      Real dN[3];
    algorithm
      (xIndex, zIndex, inArea, localX, localZ) := get_rectangle(x = x, z = z, nu = nu, nv = nv, LengthX = LengthX, LengthZ = LengthZ, OffsetX = OffsetX, OffsetZ = OffsetZ);
      if inArea and not FlatSurface then
        rectangleX := LengthX / (nu - 1);
        rectangleZ := LengthZ / (nv - 1);
        rectangleY := yMatrix[xIndex:xIndex + 1, zIndex:zIndex + 1];
        rectangleDy_dx := dy_dxMatrix[xIndex:xIndex + 1, zIndex:zIndex + 1];
        rectangleDy_dz := dy_dzMatrix[xIndex:xIndex + 1, zIndex:zIndex + 1];
        G := [rectangleY, rectangleDy_dz * rectangleZ; rectangleDy_dx * rectangleX, zeros(2, 2)];
        localXNorm := localX / rectangleX;
        localZNorm := localZ / rectangleZ;
        X := {localXNorm ^ 3, localXNorm ^ 2, localXNorm, 1};
        Z := {localZNorm ^ 3, localZNorm ^ 2, localZNorm, 1};
        dy_dx := dH * X * G * (H * Z) / rectangleX;
        dy_dz := H * X * G * (dH * Z) / rectangleZ;
        dN := cross({0, dy_dz, 1}, {1, dy_dx, 0});
        eN := dN / sqrt(dN * dN);
      else
        eN := {0, 1, 0};
      end if;
    end get_eN_protected;

    function get_elevation_protected
      extends get_elevation_base;
      input Boolean FlatSurface;
      input Modelica.Units.SI.Distance yMatrix[:, :];
      input Real dy_dxMatrix[:, :];
      input Real dy_dzMatrix[:, :];
      input Integer nu;
      input Integer nv;
      input Modelica.Units.SI.Length LengthX;
      input Modelica.Units.SI.Length LengthZ;
      input Modelica.Units.SI.Distance OffsetX;
      input Modelica.Units.SI.Distance OffsetZ;
    protected
      Integer xIndex;
      Integer zIndex;
      Modelica.Units.SI.Distance localX;
      Modelica.Units.SI.Distance localZ;
      Boolean inArea;
      Modelica.Units.SI.Distance rectangleY[2, 2];
      Real rectangleDy_dx[2, 2];
      Real rectangleDy_dz[2, 2];
      Real X[4];
      Real Z[4];
      Real G[4, 4];
      Modelica.Units.SI.Distance rectangleX;
      Modelica.Units.SI.Distance rectangleZ;
      Real localXNorm;
      Real localZNorm;
      constant Real H[4, 4] = [2, -3, 0, 1; -2, 3, 0, 0; 1, -2, 1, 0; 1, -1, 0, 0];
      //"Matrix used for the interpolation";
    algorithm
      (xIndex, zIndex, inArea, localX, localZ) := get_rectangle(x = x, z = z, nu = nu, nv = nv, LengthX = LengthX, LengthZ = LengthZ, OffsetX = OffsetX, OffsetZ = OffsetZ);
      if inArea and not FlatSurface then
        rectangleX := LengthX / (nu - 1);
        rectangleZ := LengthZ / (nv - 1);
        rectangleY := yMatrix[xIndex:xIndex + 1, zIndex:zIndex + 1];
        rectangleDy_dx := dy_dxMatrix[xIndex:xIndex + 1, zIndex:zIndex + 1];
        rectangleDy_dz := dy_dzMatrix[xIndex:xIndex + 1, zIndex:zIndex + 1];
        G := [rectangleY, rectangleDy_dz * rectangleZ; rectangleDy_dx * rectangleX, zeros(2, 2)];
        localXNorm := localX / rectangleX;
        localZNorm := localZ / rectangleZ;
        X := {localXNorm ^ 3, localXNorm ^ 2, localXNorm, 1};
        Z := {localZNorm ^ 3, localZNorm ^ 2, localZNorm, 1};
        elevation := H * X * G * (H * Z);
      else
        elevation := 0;
      end if;
    end get_elevation_protected;

    function get_surface_char_protected
      extends get_surface_char_base;
      //additional input(s)
      input Integer nu;
      input Integer nv;
      input Modelica.Units.SI.Length LengthX;
      input Modelica.Units.SI.Length LengthZ;
      input Modelica.Units.SI.Distance OffsetX;
      input Modelica.Units.SI.Distance OffsetZ;
      //"the friction coefficient";
    protected
      Integer xIndex;
      Integer zIndex;
      Boolean inArea;
    algorithm
      (xIndex, zIndex, inArea) := get_rectangle(x = x, z = z, nu = nu, nv = nv, LengthX = LengthX, LengthZ = LengthZ, OffsetX = OffsetX, OffsetZ = OffsetZ);
      if inArea and not DefaultCharSurface then
        mu := muMatrix[xIndex, zIndex];
        kc := kcMatrix[xIndex, zIndex];
        kf := kfMatrix[xIndex, zIndex];
        kv := kvMatrix[xIndex, zIndex];
        n := nMatrix[xIndex, zIndex];
        c := cMatrix[xIndex, zIndex];
        tanf := tanfMatrix[xIndex, zIndex];
        K := KMatrix[xIndex, zIndex];
        gammas := gammasMatrix[xIndex, zIndex];
      else
        mu := 0.7;
        kc := 0;
        kf := 1e10;
        kv := 1e5;
        n := 1;
        c := 2.2539e6;
        tanf := Modelica.Math.tan((44.5/180)*Modelica.Constants.pi);
        K := 0.0002;
        gammas := 14126.4;
      end if;
    end get_surface_char_protected;

    function SurfaceFunction
       extends
        Modelica.Mechanics.MultiBody.Interfaces.partialSurfaceCharacteristic;
      input Real x_min "Minimum value of x";
      input Real x_max "Maximum value of x";
      input Real z_min "Minimum value of z";
      input Real z_max "Maximum value of z";
      //input Real wz "Factor for angular frequency";
    protected
      Real aux_x;
      //Real A=(z_max-z_min)/2;
    algorithm
      for
      i in 1:nu loop
        aux_x := x_min + (x_max - x_min) * (i - 1) / (nu - 1);
        for j in 1:nv loop
          X[i, j] := aux_x;
          Z[i, j] := z_min + (z_max - z_min) * (j - 1) / (nv - 1);
          Y[i, j] := get_elevation(X[i, j], Z[i, j]);
        end for;
      end for;
      //y[i,j];
      if multiColoredSurface then
          C := {{(Y[i,j]+1)*200,255,0} for j in 1:nv, i in 1:nu};
      end if;
    end SurfaceFunction;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end Surface_base;

  model Surface_with_step
    extends Surface_base(   redeclare function get_elevation =
          get_elevation_step (                                                    FlatSurface = FlatSurface, nu = nu1, nv = nv1, LengthX = LengthX, LengthZ = LengthZ, OffsetX = OffsetX, OffsetZ = OffsetZ), redeclare
        replaceable function                                                                                                                                                                                             surfaceCharacteristic =
          SurfaceFunction_step (                                                                                                                                                                                                        x=transpose([OffsetX + 1e-6, 0.3 - 1e-3 - 1e-6, 0.3 - 1e-3 + 1e-6, 0.3 - 1e-5, 0.3 - 1e-7, 0.3 + 1e-7, 0.3 + 1e-5, 0.3 + 1e-3 - 1e-6, 0.3 + 1e-3 + 1e-6, LengthX+OffsetX - 1e-6;OffsetX + 1e-6, 0.3 - 1e-3 - 1e-6, 0.3 - 1e-3 + 1e-6, 0.3 - 1e-5, 0.3 - 1e-7, 0.3 + 1e-7, 0.3 + 1e-5, 0.3 + 1e-3 - 1e-6, 0.3 + 1e-3 + 1e-6, LengthX+OffsetX - 1e-6]), z=transpose([OffsetZ + 1e-6, OffsetZ + 1e-6, OffsetZ + 1e-6, OffsetZ + 1e-6, OffsetZ + 1e-6, OffsetZ + 1e-6, OffsetZ + 1e-6, OffsetZ + 1e-6, OffsetZ + 1e-6, OffsetZ + 1e-6;LengthZ+OffsetZ - 1e-6, LengthZ+OffsetZ - 1e-6, LengthZ+OffsetZ - 1e-6, LengthZ+OffsetZ - 1e-6, LengthZ+OffsetZ - 1e-6, LengthZ+OffsetZ - 1e-6, LengthZ+OffsetZ - 1e-6, LengthZ+OffsetZ - 1e-6, LengthZ+OffsetZ - 1e-6, LengthZ+OffsetZ - 1e-6])));
    function get_elevation_step
      extends get_elevation_base;
      input Boolean FlatSurface;
      input Integer nu;
      input Integer nv;
      input Modelica.Units.SI.Length LengthX;
      input Modelica.Units.SI.Length LengthZ;
      input Modelica.Units.SI.Distance OffsetX;
      input Modelica.Units.SI.Distance OffsetZ;
    protected
      Integer xIndex;
      Integer zIndex;
      Modelica.Units.SI.Distance localX;
      Modelica.Units.SI.Distance localZ;
      Boolean inArea;
    algorithm
      (xIndex, zIndex, inArea, localX, localZ) := get_rectangle(x = x, z = z, nu = nu, nv = nv, LengthX = LengthX, LengthZ = LengthZ, OffsetX = OffsetX, OffsetZ = OffsetZ);
      if inArea and not FlatSurface then
        elevation := (Modelica.Math.tanh(1e4*(x-0.3))+1)*0.5*0.153;
      else
        elevation := 0;
      end if;
    end get_elevation_step;

    function SurfaceFunction_step
       extends
        Modelica.Mechanics.MultiBody.Interfaces.partialSurfaceCharacteristic;

     input Real x[nu,nv] "Minimum value of x";
     input Real z[nu,nv] "Minimum value of z";
     //input Real wz "Factor for angular frequency";
    protected
      Real Y1[nu, nv];
    algorithm
      for
      i in 1:nu loop

        for j in 1:nv loop
          X[i, j] := x[i,j];
          Z[i, j] := z[i,j];
          Y[i, j] := (sign(X[i, j]-0.3)+1)*0.5*0.153;
          Y1[i, j] := ((sign(X[i, j]-(0.3-1e-3))+1)*0.5 - (sign(X[i, j]-(0.3-1e-6))+1)*0.5) + ((sign(X[i, j]-(0.3+1e-6))+1)*0.5 - (sign(X[i, j]-(0.3+1e-3))+1)*0.5);
        end for;
      end for;
      //y[i,j];
      if multiColoredSurface then
          //C := {{(Y[i,j]+1)*200,255,0} for j in 1:nv, i in 1:nu};
          C := {{0,255*(1-Y1[i,j]),0} for j in 1:nv, i in 1:nu};
      end if;
    end SurfaceFunction_step;

  end Surface_with_step;

  model contact_terrain_grouser
    final parameter Integer nx = 1;
    parameter Integer nz = 1;
    parameter Modelica.Units.SI.Length len = (0.0586*38)/80;
    parameter Modelica.Units.SI.Length lar = 0.15;
    parameter Modelica.Units.SI.Length thick = 0.012;
    parameter Real Er = 1e7;//3e6;
    parameter Real kvr = 1e3;

    parameter  Modelica.Units.SI.Position cm1 = 0;
    parameter  Modelica.Units.SI.Position cm2 = 0;
    parameter Modelica.Units.SI.Length hg = 0.0125, lg = 0;
    parameter Boolean there_is_grouser = false;
    parameter Boolean fr = false;
    parameter Boolean fl = false;
    final parameter Modelica.Units.SI.Area Ag = lg*lar / nz;
    final parameter Modelica.Units.SI.Area Af = (len * lar / (nx * nz)) - Ag;
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b annotation (
      Placement(visible = true, transformation(origin = {-100, 0}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-100, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Units.SI.Position dr2[nx * nz, 3], dr12[nx * nz, 3];
    Modelica.Units.SI.Position d2[nx * nz, 3], d12[nx * nz, 3], d3[nx * nz, 3];
    Modelica.Units.SI.Velocity vjx[nx * nz], vjz[nx * nz], modv[nx * nz], vjx1[nx * nz], vjz1[nx * nz], modv1[nx * nz], velZ[nx * nz], veldl[nx * nz], velZ1[nx * nz], veldl1[nx * nz];//, v1[nx * nz], v2[nx * nz];
    //(each stateSelect=StateSelect.prefer)
    //Modelica.Units.SI.Acceleration a1[nx * nz], a2[nx * nz];
    Modelica.Units.SI.Position jx[nx * nz], jz[nx * nz], j[nx * nz], jx1[nx * nz], jz1[nx * nz], j1[nx * nz];
  //, modj[nx * nz];
    Real cos[nx * nz], sin[nx * nz], cos1[nx * nz], sin1[nx * nz], cos2[nx * nz];
    //, cos1[nx * nz], sin1[nx * nz];
    Modelica.Units.SI.Distance deltaZ[nx * nz], deltaZ1[nx * nz], dl[nx * nz], dl1[nx * nz], dh, b[nx * nz];
    Modelica.Units.SI.Pressure p[nx * nz], preg[nx * nz];//, p2[nx * nz]
    Modelica.Units.SI.Pressure t[nx * nz];
    Modelica.Units.SI.Pressure p1[nx * nz], p1reg[nx * nz];//, p12[nx * nz]
    Modelica.Units.SI.Pressure t1[nx * nz];
    Modelica.Units.SI.Force Fbr[nx], F1[nx];
    Modelica.Units.SI.Force Fbl[nx], F2[nx];
    Modelica.Units.SI.Force Fbr1[nx], F3[nx];
    Modelica.Units.SI.Force Fbl1[nx], F4[nx];
    Modelica.Units.SI.Force Fgr[nx*nz];
    Modelica.Units.SI.Force Fgl[nx*nz];
    //, cos1[nx * nz], sin1[nx * nz];
    //Real smooth11[nx], smooth12[nx];
    Modelica.Units.SI.Force Fn, Ft, Fl;
    Modelica.Units.SI.Torque T1, T2, T3;

    Real mu[nx * nz,2];
    Real kc[nx * nz,2];
    Real kf[nx * nz,2];
    Real k[nx * nz,2];//1646666.667;
    //290000;1e10;//
    Real kv[nx * nz,2];//5e4;
    //290000;1e5;//
    Real nn[nx * nz,2];//0.13;
    ///2;1;//
    Modelica.Units.SI.ShearStress c[nx * nz,2];//700000;//0;
    //70000 * 0.5;
    Real tanf[nx * nz,2];//0.67
    Modelica.Units.SI.Length K[nx * nz,2];
    Real gammas[nx * nz,2];
    Real Nf[nx * nz,2];
    // parameter Real w = 10;

   function smooth1
      input Modelica.Units.SI.Distance d1;
      input Real w = 1e3;
      output Real s;
   algorithm
      //if d1 == 0 then
      //s := 0;
      //else
      s := (Modelica.Math.tanh(w * d1)  + 1) / 2;//- 1e-3 //(sign(vjz[nz * n])+1)*0.5 (Modelica.Math.tanh(w * (d1))  + 1) / 2
      //end if;//(Modelica.Math.tanh(w * (d1 - 1e-3)) + 1) / 2;//((Modelica.Math.atan(w * (d1 - 1e-3)) + Modelica.Constants.pi / 2) / Modelica.Constants.pi);//
      //((Modelica.Math.atan(w * d1) + Modelica.Constants.pi / 2) / Modelica.Constants.pi);
   end smooth1;

    outer AGILEX_BUNKER_PRO.Surface_base Surface;

  equation
  for n in 1:nx loop
      for m in 1:nz loop
        d2[(n - 1) * nz + m, :] = if there_is_grouser then {cm1, 0, lar * (m / nz) - lar * (0.5 / nz)} else {len * (n / nx) - len * (0.5 / nx), 0, lar * (m / nz) - lar * (0.5 / nz)};
        d12[(n - 1) * nz + m, :] = if there_is_grouser then {cm2, -hg, lar * (m / nz) - lar * (0.5 / nz)} else {0, 0, 0};
        dr2[(n - 1) * nz + m, :] = frame_b.r_0 + Modelica.Mechanics.MultiBody.Frames.resolve1(frame_b.R, d2[(n - 1) * nz + m, :]);
        dr12[(n - 1) * nz + m, :] = if there_is_grouser then frame_b.r_0 + Modelica.Mechanics.MultiBody.Frames.resolve1(frame_b.R, d12[(n - 1) * nz + m, :]) else {0, 0, 0};
        d3[(n - 1) * nz + m, :] = Modelica.Mechanics.MultiBody.Frames.resolve1(frame_b.R, {0,1,0});
        cos2[(n - 1) * nz + m] = max(1e-6, {0,1,0}*d3[(n - 1) * nz + m, :]);

        velZ[(n - 1) * nz + m] = (deltaZ[(n - 1) * nz + m] - delay(deltaZ[(n - 1) * nz + m], 1e-3))/1e-3;
        veldl[(n - 1) * nz + m] = (dl[(n - 1) * nz + m] - delay(dl[(n - 1) * nz + m], 1e-3))/1e-3;
        velZ1[(n - 1) * nz + m] = (deltaZ1[(n - 1) * nz + m] - delay(deltaZ1[(n - 1) * nz + m], 1e-3))/1e-3;
        veldl1[(n - 1) * nz + m] = (dl1[(n - 1) * nz + m] - delay(dl1[(n - 1) * nz + m], 1e-3))/1e-3;

        //a1[(n - 1) * nz + m] = (veldl[(n - 1) * nz + m] - delay(veldl[(n - 1) * nz + m], 1e-6))/1e-6;
        //a2[(n - 1) * nz + m] = (veldl1[(n - 1) * nz + m] - delay(veldl1[(n - 1) * nz + m], 1e-6))/1e-6;

        vjx[(n - 1) * nz + m] = der(dr2[(n - 1) * nz + m, :]) * Modelica.Mechanics.MultiBody.Frames.resolve1(frame_b.R, {1, 0, 0});
        vjz[(n - 1) * nz + m] = der(dr2[(n - 1) * nz + m, :]) * Modelica.Mechanics.MultiBody.Frames.resolve1(frame_b.R, {0, 0, 1});

        (mu[(n - 1) * nz + m, 1], kc[(n - 1) * nz + m, 1], kf[(n - 1) * nz + m, 1], kv[(n - 1) * nz + m, 1], nn[(n - 1) * nz + m, 1], c[(n - 1) * nz + m, 1], tanf[(n - 1) * nz + m, 1], K[(n - 1) * nz + m, 1], gammas[(n - 1) * nz + m, 1]) = Surface.get_surface_char(dr2[(n - 1) * nz + m, 1], dr2[(n - 1) * nz + m, 3]);
        (mu[(n - 1) * nz + m, 2], kc[(n - 1) * nz + m, 2], kf[(n - 1) * nz + m, 2], kv[(n - 1) * nz + m, 2], nn[(n - 1) * nz + m, 2], c[(n - 1) * nz + m, 2], tanf[(n - 1) * nz + m, 2], K[(n - 1) * nz + m, 2], gammas[(n - 1) * nz + m, 2]) = Surface.get_surface_char(dr12[(n - 1) * nz + m, 1], dr12[(n - 1) * nz + m, 3]);
        Nf[(n - 1) * nz + m,1] = Modelica.Math.tan(Modelica.Constants.pi / 4 + 0.5 * Modelica.Math.atan(tanf[(n - 1) * nz + m, 1])) ^ 2;
        Nf[(n - 1) * nz + m,2] = Modelica.Math.tan(Modelica.Constants.pi / 4 + 0.5 * Modelica.Math.atan(tanf[(n - 1) * nz + m, 2])) ^ 2;
        k[(n - 1) * nz + m, 1] = kc[(n - 1) * nz + m, 1]/(lar) + kf[(n - 1) * nz + m, 1];
        k[(n - 1) * nz + m, 2] = kc[(n - 1) * nz + m, 2]/(lar) + kf[(n - 1) * nz + m, 2];

          //v1[(n - 1) * nz + m] = -der(dr2[(n - 1) * nz + m, 2]);
          //a1[(n - 1) * nz + m] = der(v1[(n - 1) * nz + m]);
          if noEvent(dr2[(n - 1) * nz + m, 2] - Surface.get_elevation(dr2[(n - 1) * nz + m, 1], dr2[(n - 1) * nz + m, 3]) < 0) then

          der(jx[(n - 1) * nz + m]) = vjx[(n - 1) * nz + m];
          der(jz[(n - 1) * nz + m]) = vjz[(n - 1) * nz + m];
          deltaZ[(n - 1) * nz + m] + dl[(n - 1) * nz + m]*cos2[(n - 1) * nz + m] = -(dr2[(n - 1) * nz + m, 2] - Surface.get_elevation(dr2[(n - 1) * nz + m, 1], dr2[(n - 1) * nz + m, 3]));
          //p[(n - 1) * nz + m] = - smooth(0, max(0, k * (deltaZ[(n - 1) * nz + m]) ^ n - kv * der(dr2[(n - 1) * nz + m, 2])*x1));// - kv * der(dr2[(n - 1) * nz + m, 2]) smooth(0, max(0,
          //p[(n - 1) * nz + m] = - smooth(0, max(0, ((Er/(thick*0.5))*dl[(n - 1) * nz + m]) - kvr * der(dr2[(n - 1) * nz + m, 2]/cos2[(n - 1) * nz + m])*x2));// /cos2[(n - 1) * nz + m]
          (abs(preg[(n - 1) * nz + m])/k[(n - 1) * nz + m, 1])^(1/nn[(n - 1) * nz + m, 1]) = abs(deltaZ[(n - 1) * nz + m]);
          p[(n - 1) * nz + m] = - smooth(0, max(0, abs(preg[(n - 1) * nz + m])  + kv[(n - 1) * nz + m, 1] * velZ[(n - 1) * nz + m]));// - kv * der(dr2[(n - 1) * nz + m, 2]) smooth(0, max(0,
          //p2[(n - 1) * nz + m] = - smooth(0, max(0, k[(n - 1) * nz + m, 1] * abs(deltaZ[(n - 1) * nz + m]) ^ nn[(n - 1) * nz + m, 1] + kv[(n - 1) * nz + m, 1] * smooth(0, max(0, velZ[(n - 1) * nz + m]))));// - kv * der(dr2[(n - 1) * nz + m, 2]) smooth(0, max(0,
          p[(n - 1) * nz + m] = - smooth(0, max(0, ((Er/(thick*0.5))*dl[(n - 1) * nz + m]) + kvr * veldl[(n - 1) * nz + m]));// /cos2[(n - 1) * nz + m]

          //p[(n - 1) * nz + m] - p2[(n - 1) * nz + m] = 1e-3*a1[(n - 1) * nz + m];
        else

          der(jx[(n - 1) * nz + m]) = -1000 * jx[(n - 1) * nz + m];
          der(jz[(n - 1) * nz + m]) = -1000 * jz[(n - 1) * nz + m];
          deltaZ[(n - 1) * nz + m] = 0;
          dl[(n - 1) * nz + m] = 0;
          //p[(n - 1) * nz + m] - p2[(n - 1) * nz + m] = 0;
          p[(n - 1) * nz + m] = 0;
          //v1[(n - 1) * nz + m] = 0;
          //a1[(n - 1) * nz + m] = 0;
          preg[(n - 1) * nz + m] = 0;
        end if;
        j[(n - 1) * nz + m] = sqrt(jx[(n - 1) * nz + m] ^ 2 + jz[(n - 1) * nz + m] ^ 2);
        //sin1[(n - 1) * nz + m] = vjz[(n - 1) * nz + m] / (Modelica.Constants.small + sqrt(vjx[(n - 1) * nz + m] ^ 2 + vjz[(n - 1) * nz + m] ^ 2));
        //cos1[(n - 1) * nz + m] = vjx[(n - 1) * nz + m] / (Modelica.Constants.small + sqrt(vjx[(n - 1) * nz + m] ^ 2 + vjz[(n - 1) * nz + m] ^ 2));
        modv[(n - 1) * nz + m] = smooth(0, max(1e-2, sqrt(vjx[(n - 1) * nz + m] ^ 2 + vjz[(n - 1) * nz + m] ^ 2)));
        sin[(n - 1) * nz + m] = vjz[(n - 1) * nz + m] / modv[(n - 1) * nz + m];
        cos[(n - 1) * nz + m] = vjx[(n - 1) * nz + m] / modv[(n - 1) * nz + m];
        //p[(n - 1) * nz + m] = -max(0, k * deltaZ[(n - 1) * nz + m] ^ n + kv * der(deltaZ[(n - 1) * nz + m]));
        t[(n - 1) * nz + m] = min(-mu[(n - 1) * nz + m, 1] * p[(n - 1) * nz + m],(c[(n - 1) * nz + m, 1] - tanf[(n - 1) * nz + m, 1] * p[(n - 1) * nz + m])) * (1 - exp(-j[(n - 1) * nz + m] / K[(n - 1) * nz + m, 1]));
        //modj[(n - 1) * nz + m] = max(1e-6, j[(n - 1) * nz + m]);
        //sin[(n - 1) * nz + m] = jz[(n - 1) * nz + m] / modj[(n - 1) * nz + m];
        //cos[(n - 1) * nz + m] = jx[(n - 1) * nz + m] / modj[(n - 1) * nz + m];
        //cos[(n - 1) * nz + m] = (Modelica.Math.tanh(1e3 * vjx[(n - 1) * nz + m])) * sqrt(1 - sin[(n - 1) * nz + m]^2);

        vjx1[(n - 1) * nz + m] = der(dr12[(n - 1) * nz + m, :]) * Modelica.Mechanics.MultiBody.Frames.resolve1(frame_b.R, {1, 0, 0});
        vjz1[(n - 1) * nz + m] = der(dr12[(n - 1) * nz + m, :]) * Modelica.Mechanics.MultiBody.Frames.resolve1(frame_b.R, {0, 0, 1});
        // v2[(n - 1) * nz + m] = -der(dr12[(n - 1) * nz + m, 2]);
        //a2[(n - 1) * nz + m] = der(v2[(n - 1) * nz + m]);
        if noEvent(there_is_grouser and dr12[(n - 1) * nz + m, 2] - Surface.get_elevation(dr12[(n - 1) * nz + m, 1], dr12[(n - 1) * nz + m, 3]) < 0) then

          der(jx1[(n - 1) * nz + m]) = vjx1[(n - 1) * nz + m];
          der(jz1[(n - 1) * nz + m]) = vjz1[(n - 1) * nz + m];

          //p1[(n - 1) * nz + m] = - smooth(0, max(0, k * (deltaZ1[(n - 1) * nz + m]) ^ n - kv * der(dr12[(n - 1) * nz + m, 2])*x1));// - kv * der(dr12[(n - 1) * nz + m, 2]) smooth(0, max(0,

        //deltaZ1[(n - 1) * nz + m] + dl1[(n - 1) * nz + m]*cos2[(n - 1) * nz + m] = - (dr12[(n - 1) * nz + m, 2] - Surface.get_elevation(dr12[(n - 1) * nz + m, 1], dr12[(n - 1) * nz + m, 3]));

          //p1[(n - 1) * nz + m] = - smooth(0, max(0, ((Er/(thick*0.5 + hg))*dl1[(n - 1) * nz + m]) - kvr * der(dr12[(n - 1) * nz + m, 2]/cos2[(n - 1) * nz + m])*x2));// /cos2[(n - 1) * nz + m]

          (abs(p1reg[(n - 1) * nz + m])/ k[(n - 1) * nz + m, 2])^(1/nn[(n - 1) * nz + m, 2]) = abs(deltaZ1[(n - 1) * nz + m]);
          p1[(n - 1) * nz + m] = - smooth(0, max(0, abs(p1reg[(n - 1) * nz + m]) + kv[(n - 1) * nz + m, 2] * velZ1[(n - 1) * nz + m]));// - kv * der(dr12[(n - 1) * nz + m, 2]) smooth(0, max(0,
          //p12[(n - 1) * nz + m] = - smooth(0, max(0, k[(n - 1) * nz + m, 2] * abs(deltaZ1[(n - 1) * nz + m]) ^ nn[(n - 1) * nz + m, 2] + kv[(n - 1) * nz + m, 2] * smooth(0, max(0, velZ1[(n - 1) * nz + m]))));// - kv * der(dr12[(n - 1) * nz + m, 2]) smooth(0, max(0,
        deltaZ1[(n - 1) * nz + m] + dl1[(n - 1) * nz + m]*cos2[(n - 1) * nz + m] = - (dr12[(n - 1) * nz + m, 2] - Surface.get_elevation(dr12[(n - 1) * nz + m, 1], dr12[(n - 1) * nz + m, 3]));

          p1[(n - 1) * nz + m] = - smooth(0, max(0, ((Er/(thick*0.5 + hg))*dl1[(n - 1) * nz + m]) + kvr * veldl1[(n - 1) * nz + m]));// /cos2[(n - 1) * nz + m]

        //p1[(n - 1) * nz + m] - p12[(n - 1) * nz + m] = 1e-3*a2[(n - 1) * nz + m];
        else

          der(jx1[(n - 1) * nz + m]) = -1000 * jx1[(n - 1) * nz + m];
          der(jz1[(n - 1) * nz + m]) = -1000 * jz1[(n - 1) * nz + m];
          deltaZ1[(n - 1) * nz + m] = 0;
          dl1[(n - 1) * nz + m] = 0;
          //p1[(n - 1) * nz + m] - p12[(n - 1) * nz + m] = 0;
          p1[(n - 1) * nz + m] = 0;
          //v2[(n - 1) * nz + m] = 0;
          //a2[(n - 1) * nz + m] = 0;
          p1reg[(n - 1) * nz + m] = 0;
        end if;
        j1[(n - 1) * nz + m] = sqrt(jx1[(n - 1) * nz + m] ^ 2 + jz1[(n - 1) * nz + m] ^ 2);
        //sin1[(n - 1) * nz + m] = vjz[(n - 1) * nz + m] / (Modelica.Constants.small + sqrt(vjx[(n - 1) * nz + m] ^ 2 + vjz[(n - 1) * nz + m] ^ 2));
        //cos1[(n - 1) * nz + m] = vjx[(n - 1) * nz + m] / (Modelica.Constants.small + sqrt(vjx[(n - 1) * nz + m] ^ 2 + vjz[(n - 1) * nz + m] ^ 2));
        modv1[(n - 1) * nz + m] = smooth(0, max(1e-2, sqrt(vjx1[(n - 1) * nz + m] ^ 2 + vjz1[(n - 1) * nz + m] ^ 2)));
        sin1[(n - 1) * nz + m] = vjz1[(n - 1) * nz + m] / modv1[(n - 1) * nz + m];
        cos1[(n - 1) * nz + m] = vjx1[(n - 1) * nz + m] / modv1[(n - 1) * nz + m];
        //p[(n - 1) * nz + m] = -max(0, k * deltaZ[(n - 1) * nz + m] ^ n + kv * der(deltaZ[(n - 1) * nz + m]));
        t1[(n - 1) * nz + m] = min(-mu[(n - 1) * nz + m, 2] * p1[(n - 1) * nz + m],(c[(n - 1) * nz + m, 2] - tanf[(n - 1) * nz + m, 2] * p1[(n - 1) * nz + m])) * (1 - exp(-j1[(n - 1) * nz + m] / K[(n - 1) * nz + m, 2]));
        //modj[(n - 1) * nz + m] = max(1e-6, j[(n - 1) * nz + m]);
        //sin[(n - 1) * nz + m] = jz[(n - 1) * nz + m] / modj[(n - 1) * nz + m];
        //cos[(n - 1) * nz + m] = jx[(n - 1) * nz + m] / modj[(n - 1) * nz + m];
        //cos[(n - 1) * nz + m] = (Modelica.Math.tanh(1e3 * vjx[(n - 1) * nz + m])) * sqrt(1 - sin[(n - 1) * nz + m]^2);
        Fgl[(n - 1) * nz + m] = if fl and noEvent((deltaZ1[(n - 1) * nz + m]/cos2[(n - 1) * nz + m]) > 1e-3) then -(0.5 * gammas[(n - 1) * nz + m, 2] * min(deltaZ1[(n - 1) * nz + m]/cos2[(n - 1) * nz + m], hg) ^ 2 * Nf[(n - 1) * nz + m, 2] - p[(n - 1) * nz + m]*min(deltaZ1[(n - 1) * nz + m]/cos2[(n - 1) * nz + m], hg)*Nf[(n - 1) * nz + m, 2] + 2 * c[(n - 1) * nz + m, 2] * min(deltaZ1[(n - 1) * nz + m]/cos2[(n - 1) * nz + m], hg) * sqrt(Nf[(n - 1) * nz + m, 2])) * (lar / nz) * (1 - exp(-abs(jx1[(n - 1) * nz + m]) / K[(n - 1) * nz + m, 2])) * smooth1(-vjx1[(n - 1) * nz + m]) else 0; //(sign(-vjx1[(n - 1) * nz + m])+1)*0.5
        Fgr[(n - 1) * nz + m] = if fr and noEvent((deltaZ1[(n - 1) * nz + m]/cos2[(n - 1) * nz + m]) > 1e-3) then (0.5 * gammas[(n - 1) * nz + m, 2] * min(deltaZ1[(n - 1) * nz + m]/cos2[(n - 1) * nz + m], hg) ^ 2 * Nf[(n - 1) * nz + m, 2] - p[(n - 1) * nz + m]*min(deltaZ1[(n - 1) * nz + m]/cos2[(n - 1) * nz + m], hg)*Nf[(n - 1) * nz + m, 2] + 2 * c[(n - 1) * nz + m, 2] * min(deltaZ1[(n - 1) * nz + m]/cos2[(n - 1) * nz + m], hg) * sqrt(Nf[(n - 1) * nz + m, 2])) * (lar / nz) * (1 - exp(-abs(jx1[(n - 1) * nz + m]) / K[(n - 1) * nz + m, 2])) * smooth1(vjx1[(n - 1) * nz + m]) else 0; //(sign(vjx1[(n - 1) * nz + m])+1)*0.5
        b[(n - 1) * nz + m] = (-d12[(n - 1) * nz + m, 2] - dl1[(n - 1) * nz + m] - min(deltaZ1[(n - 1) * nz + m]/cos2[(n - 1) * nz + m], hg  - dl1[(n - 1) * nz + m]) *0.5);
      end for;
      F1[n] = if noEvent((deltaZ[nz * n]/cos2[nz * n]) > 1e-3) then (0.5 * gammas[nz * n, 1] * min(deltaZ[nz * n]/cos2[nz * n], thick) ^ 2 * Nf[nz * n, 1] + 2 * c[nz * n, 1] * min(deltaZ[nz * n]/cos2[nz * n], thick) * sqrt(Nf[nz * n, 1])) * (len / nx) else 0;
      F2[n] = if noEvent((deltaZ[nz * (n - 1) + 1]/cos2[nz * (n - 1) + 1]) > 1e-3) then (0.5 * gammas[nz * (n - 1) + 1, 1] * min(deltaZ[nz * (n - 1) + 1]/cos2[nz * (n - 1) + 1], thick) ^ 2 * Nf[nz * (n - 1) + 1, 1] + 2 * c[nz * (n - 1) + 1, 1] * min(deltaZ[nz * (n - 1) + 1]/cos2[nz * (n - 1) + 1], thick) * sqrt(Nf[nz * (n - 1) + 1, 1])) * (len / nx) else 0;
      F3[n] = if noEvent((deltaZ1[nz * n]/cos2[nz * n]) > 1e-3) then (0.5 * gammas[nz * n, 2] * min(deltaZ1[nz * n]/cos2[nz * n],hg) ^ 2 * Nf[nz * n, 2] + 2 * c[nz * n, 2] * min(deltaZ1[nz * n]/cos2[nz * n],hg) * sqrt(Nf[nz * n, 2])) * (lg) else 0;
      F4[n] = if noEvent((deltaZ1[nz * (n - 1) + 1]/cos2[nz * (n - 1) + 1]) > 1e-3) then (0.5 * gammas[nz * (n - 1) + 1, 2] * min(deltaZ1[nz * (n - 1) + 1]/cos2[nz * (n - 1) + 1],hg) ^ 2 * Nf[nz * (n - 1) + 1, 2] + 2 * c[nz * (n - 1) + 1, 2] *min(deltaZ1[nz * (n - 1) + 1]/cos2[nz * (n - 1) + 1],hg) * sqrt(Nf[nz * (n - 1) + 1, 2])) * (lg) else 0;
      Fbr[n] = min((abs(vjz[nz * n])/1e-2)*F1[n], F1[n]) * smooth1(vjz[nz * n]);//smooth1(vjz[nz * n]); (sign(vjz[nz * n])+1)*0.5
      Fbl[n] = min((abs(vjz[nz * (n - 1) + 1])/1e-2)*F2[n], F2[n]) * (-smooth1(-vjz[nz * (n - 1) + 1]));//(-smooth1(-vjz[nz * (n - 1) + 1])); (-(sign(-vjz[nz * (n - 1) + 1])+1)*0.5)
      Fbr1[n] = min((abs(vjz1[nz * n])/1e-2)*F3[n], F3[n]) * smooth1(vjz1[nz * n]);//smooth1(vjz1[nz * n]); (sign(vjz1[nz * n])+1)*0.5
      Fbl1[n] = min((abs(vjz1[nz * (n - 1) + 1])/1e-2)*F4[n], F4[n]) * (-smooth1(-vjz1[nz * (n - 1) + 1]));//(-smooth1(-vjz1[nz * (n - 1) + 1])); (-(sign(-vjz1[nz * (n - 1) + 1])+1)*0.5)
    end for;
    Ft = t * cos * Af + t1 * cos1 * Ag  + sum(Fgr) + sum(Fgl);
    Fn = sum(p) * Af + sum(p1) * Ag;
    Fl = t * sin * Af + sum(Fbr) + sum(Fbl) + t1 * sin1 * Ag + sum(Fbr1) + sum(Fbl1);
    frame_b.f = {Ft, Fn, Fl};
    //frame_b.f = {t * cos * (len * lar / (nx * nz)), sum(p) * (len * lar / (nx * nz)), t * sin * (len * lar / (nx * nz)) + sum(Fbr) + sum(Fbl)};
    T1 = -p * d2[:, 3] * Af-p1 * d12[:, 3] * Ag - sum(t1[i] * Ag * (-d12[i, 2] * sin1[i]) for i in 1:nx * nz) - sum(Fbr1[i] * b[nz * i] for i in 1:nx) - sum(Fbl1[i] * b[nz * (i - 1) + 1] for i in 1:nx)
    + sum(Fbr[i] * min(deltaZ[nz * i]/cos2[nz * i], thick)*0.5 for i in 1:nx) + sum(Fbl[i] * min(deltaZ[nz * (i - 1) + 1]/cos2[nz * (i - 1) + 1], thick)*0.5 for i in 1:nx);
    T2 = sum(t[i] * Af * (d2[i, 3] * cos[i] - d2[i, 1] * sin[i]) for i in 1:nx * nz) - sum(Fbr[i] * d2[nz * i, 1] for i in 1:nx) - sum(Fbl[i] * d2[nz * (i - 1) + 1, 1] for i in 1:nx) +
    sum(t1[i] * Ag * (d12[i, 3] * cos1[i] - d12[i, 1] * sin1[i]) for i in 1:nx * nz) - sum(Fbr1[i] * d12[nz * i, 1] for i in 1:nx) - sum(Fbl1[i] * d12[nz * (i - 1) + 1, 1] for i in 1:nx) +
    sum((Fgl[i] + Fgr[i]) * d12[i, 3]  for i in 1:nx * nz);
    T3 = p * d2[:, 1] * Af + p1 * d12[:, 1] * Ag + sum(t1[i] * Ag * (-d12[i, 2] * cos1[i]) for i in 1:nx * nz)+ sum((Fgl[i] + Fgr[i]) * b[i]  for i in 1:nx * nz);
    frame_b.t = {T1, T2, T3};

  algorithm

    dh := sum(dl1)/(nx*nz);
    //frame_b.t = {-p * d[:, 3] * (len * lar / (nx * nz)), sum(t[i] * (len * lar / (nx * nz)) * (d[i, 3] * cos[i] - d[i, 1] * sin[i]) for i in 1:nx * nz) + sum(Fbr[i] * d[nz * i, 1] for i in 1:nx) + sum(Fbl[i] * d[nz * (i - 1) + 1, 1] for i in 1:nx), p * d[:, 1] * (len * lar / (nx * nz))};
    //for i in 1:nx loop
    //smooth11[i] = smooth1(vjz[nz * i]);
    //smooth12[i] = -smooth1(-vjz[nz * (i - 1) + 1]);
    //end for;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end contact_terrain_grouser;

  model contact_t_grouser
    parameter Integer nz = 1;
    parameter Integer n = 80;
    parameter Boolean there_is_grouser[n] = {true,true,true,true,true,false,true,false,true,false,true,false,true,false,true,
    false,true,true,true,true,true,true,true,true,false,true,false,true,false,true,false,true,false,true,false,true,true,true,
    true,true,true,true,true,true,true,false,true,false,true,false,true,false,true,false,true,false,true,true,true,true,true,
    true,true,true,false,true,false,true,false,true,false,true,false,true,false,true,true,true,true,true};
     parameter Modelica.Units.SI.Position cm1[n] = {0.0175675,0.0117325,0.0190325,0.0131975,0.0204975,0,0.019573828,0,0.01634164,0,
     0.013109453,0,0.009877265,0,0.0066450783,0,0.00807,0.01537,0.009535,0.016835,0.011,0.0183,0.012465,0.019765,0,0.021189922,0,
     0.017957734,0,0.014725547,0,0.01149336,0,0.008261172,0,0.0073375,0.0146375,0.0088025,0.0161025,0.0102675,0.0175675,0.0117325,
     0.0190325,0.0131975,0.0204975,0,0.019573828,0,0.01634164,0,0.013109453,0,0.009877265,0,0.0066450783,0,0.00807,0.01537,0.009535,
     0.016835,0.011,0.0183,0.012465,0.019765,0,0.021189922,0,0.017957734,0,0.014725547,0,0.01149336,0,0.008261172,0,0.0073375,0.0146375,
     0.0088025,0.0161025,0.0102675};
      parameter Modelica.Units.SI.Position cm2[n] = {0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,0.01172,0,0.01465,0,0.01758,
      0,0.02051,0,0.0219875,0.0014525,0.0234525,0.0029175,0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,0.010255,0,0.013185,0,
      0.016115,0,0.019045,0,0.021255,0.00072,0.02272,0.002185,0.024185,0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,0.01172,0,
      0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,0.0029175,0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,0.010255,
      0,0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,0.002185,0.024185};
        parameter Modelica.Units.SI.Length lg[n] = {0.0073,0.00437,0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,
        0.011695,0.002905,0.008765,0.005835,0.005835,0.008765,0.002905,0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.01316,0.00144,
        0.01023,0.00437,0.0073,0.0073,0.00437,0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.011695,0.002905,0.008765,
        0.005835,0.005835,0.008765,0.002905,0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.01316,0.00144,0.01023,0.00437,0.0073};
       parameter Boolean fr[n]= {true,false,true,false,true,false,true,false,true,false,true,false,true,false,true,false,false,true,
       false,true,false,true,false,true,false,true,false,true,false,true,false,true,false,true,false,false,true,false,true,false,
       true,false,true,false,true,false,true,false,true,false,true,false,true,false,true,false,false,true,false,true,false,true,false,
       true,false,true,false,true,false,true,false,true,false,true,false,false,true,false,true,false};
        parameter Boolean fl[n]= {false,true,false,true,false,false,true,false,true,false,true,false,true,false,true,false,true,false,
        true,false,true,false,true,false,false,true,false,true,false,true,false,true,false,true,false,true,false,true,false,true,false,
        true,false,true,false,false,true,false,true,false,true,false,true,false,true,false,true,false,true,false,true,false,true,false,
        false,true,false,true,false,true,false,true,false,true,false,true,false,true,false,true};
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b[n] annotation (
      Placement(visible = true, transformation(origin = {-118, 3}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-118, 4}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    AGILEX_BUNKER_PRO.contact_terrain_grouser contact_terrain[n](
      there_is_grouser=there_is_grouser,
      cm1=cm1,
      cm2=cm2,
      lg=lg,
      fr=fr,
      fl=fl,
      each nz=nz) annotation (Placement(visible=true, transformation(
          origin={-2,0},
          extent={{-10,-10},{10,10}},
          rotation=0)));

  equation
    for i in 1:n loop
      connect(contact_terrain[i].frame_b, frame_b[i]);
    end for;
    annotation (
      Icon(coordinateSystem(grid = {2, 0})),
      Diagram(coordinateSystem(grid = {2, 3})));
  end contact_t_grouser;

  model script_grouser
    parameter Real larg = 0.0146;
    parameter Real distg = 0.0586;
    parameter Integer ng = 38;
    parameter Integer n = 80;
    parameter Real lar = 0.15;
    final parameter Real xg[ng] = {distg*(i-1) for i in 1:ng};

    Boolean there_is_gouser[n](each fixed = true, each start = false);
    Real cm1[n](each fixed = true, each start = 0);
    Real cm2[n](each fixed = true, each start = 0);
    Real lg[n](each fixed = true, each start = 0);
    Boolean fr[n](each fixed = true, each start = false);
    Boolean fl[n](each fixed = true, each start = false);

  protected
    Integer rectX;
    Integer rectZ;
    Boolean inArea;
    Modelica.Units.SI.Distance localX;
    Modelica.Units.SI.Distance localZ;
    Integer rectX1;
    Integer rectZ1;
    Boolean inArea1;
    Modelica.Units.SI.Distance localX1;
    Modelica.Units.SI.Distance localZ1;

  function get_rectangle
      //"finds the rectangle which calculates the elevation and normal vector for the current position of the tire."
      input Modelica.Units.SI.Position x;
      //"acutal position to find active recthangle and relative position in x direction";
      input Modelica.Units.SI.Position z;
      //"acutal position to find active recthangle and relative position in z direction";
      input Integer nu;
      // "Number of grid points in x direction";
      input Integer nv;
      // "Number of grid points in z direction";
      input Modelica.Units.SI.Length LengthX;
      // "Length of the surface area in x direction";
      input Modelica.Units.SI.Length LengthZ;
      // "Length of the surface area in z direction";
      input Modelica.Units.SI.Distance OffsetX;
      // "Offset of the surface area in x direction";
      input Modelica.Units.SI.Distance OffsetZ;
      // "Offset of the surface area in z direction";
      output Integer rectX;
      // "index of the active rectange in x direction";
      output Integer rectZ;
      // "index of the active rectange in z direction";
      output Boolean inArea;
      //"true: tire is in the defined surface area, false: it is not.";
      output Modelica.Units.SI.Distance localX;
      //"relative position in the active rectangle in x direction.";
      output Modelica.Units.SI.Distance localZ;
      //"relative position in the active rectangle in z direction.";
    protected
      Modelica.Units.SI.Distance rectangleX;
      Modelica.Units.SI.Distance rectangleZ;
  algorithm
      rectX := integer(ceil((x - OffsetX) / LengthX * (nu - 1)));
      rectZ := integer(ceil((z - OffsetZ) / LengthZ * (nv - 1)));
      if rectX > nu - 1 or rectX < 1 or rectZ > nv - 1 or rectZ < 1 then
        inArea := false;
      else
        inArea := true;
      end if;
      if rectX > nu - 1 then
        rectX :=1;
      elseif rectX < 1 then
        rectX :=nu - 1;
      end if;
      rectangleX := LengthX / (nu - 1);
      rectangleZ := LengthZ / (nv - 1);
      // getting local position (in the rectangle the tire is in at the time of calling)
      localX := mod(x - LengthX - OffsetX, rectangleX);
      localZ := mod(z - LengthZ - OffsetZ, rectangleZ);
  end get_rectangle;

  algorithm

  for i in 1:ng loop
    (rectX, rectZ, inArea, localX, localZ) := get_rectangle(xg[i] - larg/2, lar*0.5, n +1, 2, distg*ng, lar, 0, 0);
    (rectX1, rectZ1, inArea1, localX1, localZ1) := get_rectangle(xg[i] + larg/2, lar*0.5, n+1, 2, distg*ng, lar, 0, 0);

    if rectX == rectX1 then
      there_is_gouser[rectX] := true;
      cm1[rectX] := (localX*localX*0.5 + (distg*ng/n - localX1)*(distg*ng/n + localX1)*0.5)/(localX + (distg*ng/n - localX1));
      cm2[rectX] := (localX + localX1)*0.5;
      lg[rectX] := larg;
      fr[rectX] := true;
      fl[rectX] := true;
    elseif rectX1 == rectX + 1 or (rectX == n and rectX1 == 1) then
      there_is_gouser[rectX] := true;
      there_is_gouser[rectX1] := true;
      cm1[rectX] := localX*0.5;
      cm2[rectX] := (localX + distg*ng/n)*0.5;
      lg[rectX] := -localX + distg*ng/n;
      fl[rectX] := true;
      cm1[rectX1] := (localX1 + distg*ng/n)*0.5;
      cm2[rectX1] := localX1*0.5;
      lg[rectX1] := localX1;
      fr[rectX1] := true;
    end if;

  end for;

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=20,
        __Dymola_NumberOfIntervals=2500,
        Tolerance=0.001,
        __Dymola_Algorithm="Cvode"));
  end script_grouser;

  model grouser
    parameter Boolean anim = true;
    input Modelica.Units.SI.Length hg = 0.0125;
    parameter Modelica.Units.SI.Length lg = 0.0144728007;
    parameter Modelica.Units.SI.Length lar = 0.15;
    outer Modelica.Mechanics.MultiBody.World world;
    Modelica.Mechanics.MultiBody.Visualizers.Advanced.Shape shape(
      shapeType="box",
      R=frame_a.R,
      r=frame_a.r_0,
      lengthDirection={0,-1,0},
      widthDirection={0,0,1},
      length=hg,
      width=lar,
      height=lg,
      color={0,0,0},
      specularCoefficient=world.defaultSpecularCoefficient) if anim
      annotation (Placement(transformation(extent={{-6,-2},{14,18}})));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a
      annotation (Placement(transformation(extent={{-112,-10},{-80,22}})));
    Modelica.Mechanics.MultiBody.Parts.Body body(
      animation=false,
      r_CM={0,-0.0125*0.5,0},
      m=0.03285,
      I_11=(1/12)*0.03285*(0.0125^2 + 0.15^2),
      I_22=(1/12)*0.03285*(lg^2 + 0.15^2),
      I_33=(1/12)*0.03285*(lg^2 + 0.0125^2))
             annotation (Placement(transformation(extent={{-58,-6},{-38,14}})));
  equation
    connect(frame_a, body.frame_a) annotation (Line(
        points={{-96,6},{-94,6},{-94,4},{-58,4}},
        color={95,95,95},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end grouser;

  model element_grouser
    parameter Boolean anim = true;
    parameter Modelica.Units.SI.Position start[3] = {0,0,0};
  parameter Modelica.Units.SI.Position startX = 0;
    parameter Modelica.Units.SI.Position startY = 0;
    parameter Modelica.Units.SI.Position startZ = 0;
    parameter Modelica.Units.SI.Angle thZ = 0;
    parameter Modelica.Units.SI.Angle a = 0;
    parameter Modelica.Units.SI.Angle b = 0;
    parameter Modelica.Units.SI.Angle g = 0;
    final parameter Modelica.Mechanics.MultiBody.Frames.Orientation R = Modelica.Mechanics.MultiBody.Frames.axesRotations({1,2,3}, {a,b,g}, zeros(3));
    final parameter Modelica.Units.SI.Position starte[3] = start + Modelica.Mechanics.MultiBody.Frames.resolve1(R, {startX, startY, startZ});
    parameter Modelica.Units.SI.Position cm = 0;
    input Modelica.Units.SI.Length dh = 0;
    Modelica.Units.SI.Length hg = 0.0125 - dh;
    parameter Modelica.Units.SI.Length lg = 0;
    parameter Boolean there_is_grouser = false;
    parameter Real l = 2.28;
    Modelica.Mechanics.MultiBody.Parts.BodyShape bodyShape(
      animation=anim,
      I_11=(1/12)*6.2517*(0.012^2 + 0.15^2),
      I_22=(1/12)*6.2517*((l/80)^2 + 0.15^2),
      I_33=(1/12)*6.2517*((l/80)^2 + 0.012^2),             angles_fixed = true, angles_start(each displayUnit = "rad") = {a, b, thZ+g}, animateSphere = false, height = 0.15,
      length=l/80,
      m=6.2517/80,
      r={l/80,0,0},
      r_0(fixed=true, start={starte[1],starte[2],starte[3]}),
      r_CM={(l/80)*0.5,0,0},                                                                                                                                                                                                        shapeType = "box",
      color={0,0,0},                                                                                                                                                                                                        useQuaternions = false,
      width=0.012)                                                                                                                                                                                                         annotation (
      Placement(visible = true, transformation(origin = {-12, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a annotation (
      Placement(visible = true, transformation(origin = {-92, 0}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-92, 0}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a1 annotation (
      Placement(visible = true, transformation(origin = {-8, -74}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-4, -98}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation(animation = false, r={0,-0.006,
          -0.15*0.5})                                                                                                      annotation (
      Placement(visible = true, transformation(origin = {-10, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation1(animation = false, r={(l/80)
          *0.5,0.006,0})                                                                                                   annotation (
      Placement(visible = true, transformation(origin = {-26, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    AGILEX_BUNKER_PRO.bush_element bush_element annotation (Placement(visible=
            true, transformation(
          origin={38,4},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b1 annotation (
      Placement(visible = true, transformation(origin = {-12, 72}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(extent = {{-20, 76}, {12, 108}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b annotation (
      Placement(visible = true, transformation(origin = {90, 6}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Units.SI.Velocity v0[3];
    Modelica.Units.SI.Position r0[3];
    Modelica.Units.SI.Angle th0;
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation2(
        animation=false, r={cm,-0.012*0.5,0})
                                     if there_is_grouser
      annotation (Placement(transformation(extent={{-68,-64},{-48,-44}})));
    grouser gouser1(
      anim=anim,    hg=hg, lg=lg) if there_is_grouser
      annotation (Placement(transformation(extent={{-36,-66},{-16,-46}})));
  equation
    v0 = der(bodyShape.frame_a.r_0);
    r0 = bodyShape.frame_a.r_0;
    th0 = Modelica.Math.atan2(bodyShape.frame_b.r_0[2] - bodyShape.frame_a.r_0[2],bodyShape.frame_b.r_0[1] - bodyShape.frame_a.r_0[1]);
    connect(frame_a, bodyShape.frame_a) annotation (
      Line(points = {{-92, 0}, {-22, 0}, {-22, 4}}));
    connect(bodyShape.frame_a, fixedTranslation.frame_a) annotation (
      Line(points = {{-22, 4}, {-36, 4}, {-36, -32}, {-20, -32}}, color = {95, 95, 95}));
    connect(fixedTranslation1.frame_a, bodyShape.frame_a) annotation (
      Line(points = {{-36, 40}, {-56, 40}, {-56, 4}, {-22, 4}}, color = {95, 95, 95}));
    connect(bodyShape.frame_b, bush_element.frame_a) annotation (
      Line(points={{-2,4},{14,4},{14,3.8},{28.2,3.8}}));
    connect(fixedTranslation1.frame_b, frame_b1) annotation (
      Line(points = {{-16, 40}, {-12, 40}, {-12, 72}}));
    connect(bush_element.frame_b, frame_b) annotation (
      Line(points={{48,3.8},{90,3.8},{90,6}},    color = {95, 95, 95}));
    connect(fixedTranslation.frame_b, frame_a1) annotation (
      Line(points = {{0, -32}, {76, -32}, {76, -74}, {-8, -74}}));
    connect(bodyShape.frame_a, fixedTranslation2.frame_a) annotation (Line(
        points={{-22,4},{-52,4},{-52,-24},{-96,-24},{-96,-54},{-68,-54}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation2.frame_b, gouser1.frame_a) annotation (Line(
        points={{-48,-54},{-48,-55.4},{-35.6,-55.4}},
        color={95,95,95},
        thickness=0.5));
    annotation (
      Icon(graphics={  Rectangle(origin = {-6, -1}, fillColor = {113, 113, 113}, fillPattern = FillPattern.Solid, extent = {{-34, 99}, {34, -99}})}));
  end element_grouser;

  model track_with_grouser
    parameter Boolean anim = true;
    parameter Integer n = 80;
    parameter Real start[3] = {0, 0, 0};
    input Modelica.Units.SI.Length dh[n];
    parameter Real l = 0.0586*38;
    parameter Modelica.Units.SI.Length startX[n] = {0.13524441,0.16336884,0.19149324,0.21961766,0.24774207,0.27586648,0.3039909,
    0.33211532,0.3602397,0.38836285,0.41647583,0.44458714,0.47269815,0.5008085,0.52891237,0.55688405,0.58262175,0.6032808,0.62217665,
    0.6408851,0.65957004,0.6782187,0.696536,0.7112605,0.716824,0.71197385,0.69739616,0.67507076,0.64832556,0.62066644,0.59294254,0.56521285,
    0.537483,0.5097535,0.48202458,0.4542961,0.42656812,0.39884064,0.37111363,0.3433871,0.31566107,0.28793547,0.26020995,0.23248067,
    0.20471112,0.17665297,0.1495237,0.123243846,0.09706277,0.07089086,0.04471898,0.018546196,-0.007627586,-0.03380237,-0.059978165,
    -0.08615496,-0.11233275,-0.13851139,-0.1646895,-0.19085327,-0.21686272,-0.24097931,-0.25945207,-0.27013397,-0.2719339,-0.26467142,
    -0.24910323,-0.226839,-0.20030162,-0.17248595,-0.14458823,-0.11668263,-0.0887759,-0.06086873,-0.032961085,-0.0050521013,0.022866843,
    0.050872553,0.07899562,0.10712001} + fill(0.18554, n);
    parameter Modelica.Units.SI.Length startY[n] = {-0.19040328,-0.19041218,-0.19041839,-0.19042167,-0.19042195,-0.19041924,-0.19041334,
    -0.19040224,-0.19036625,-0.19010034,-0.18929689,-0.18843658,-0.1875661,-0.18667395,-0.18559812,-0.18267798,-0.17135578,-0.15227471,
    -0.1314442,-0.11044516,-0.08942516,-0.068372965,-0.04703182,-0.023075268,0.004486749,0.032183215,0.05622712,0.073321044,0.08200841,
    0.08709987,0.09182613,0.096518226,0.101209484,0.10590327,0.11059991,0.11529946,0.120001905,0.12470725,0.1294155,0.13412665,0.13884066,
    0.14355734,0.1482744,0.1529695,0.15742008,0.15932229,0.15191655,0.14190124,0.13163047,0.12133636,0.11104223,0.10075044,0.0904612,0.08017456,
    0.06989051,0.05960904,0.049330126,0.03905341,0.02877536,0.018460942,0.007763456,-0.006696874,-0.027895795,-0.05390605,-0.08196663,-0.10913077,
    -0.1325458,-0.14971939,-0.15901689,-0.16316524,-0.16672072,-0.17021371,-0.17369777,-0.17717819,-0.18065481,-0.1841206,-0.18750517,-0.19007342,
    -0.19034702,-0.19038792} + fill(0.07342, n);
    parameter Modelica.Units.SI.Length startZ[n] = fill(0, n);
     parameter Modelica.Units.SI.Angle thZ[n] = {-0.0003167978,-0.00022154418,-0.00011671506,-1.08526265e-05,9.582425e-05,0.00020965889,0.00039465772,
     0.0012881175,0.009234589,0.028549504,0.03059104,0.030956108,0.0317336,0.038335502,0.10310366,0.4122258,0.7446929,0.8339565,0.8430107,0.8441416,
     0.8458825,0.861678,1.0180945,1.3697258,1.7422643,2.1139288,2.4861243,2.8260055,2.9593987,2.972725,2.9739714,2.9740026,2.9739118,2.9738085,2.973704,
     2.9735994,2.973495,2.9733903,2.973286,2.9731822,2.9730868,2.9730735,2.973875,2.9827793,3.0722399,-2.8762295,-2.7775958,-2.7677548,-2.766853,-2.7668512,
     -2.7669399,-2.7670376,-2.7671363,-2.767235,-2.7673337,-2.7674313,-2.7675147,-2.7674623,-2.766057,-2.7512207,-2.6029167,-2.2892492,-1.9621332,-1.6365067,
     -1.3112022,-0.9856874,-0.65874875,-0.33877757,-0.14829077,-0.12678835,-0.12452726,-0.124204345,-0.12407412,-0.1239377,-0.12354662,-0.120608404,
     -0.092386425,-0.00982342,-0.0014643206,-0.00054822356};
          parameter Boolean there_is_grouser[n] = {true,true,true,true,true,false,true,false,true,false,true,false,true,false,true,false,true,true,true,
          true,true,true,true,true,false,true,false,true,false,true,false,true,false,true,false,true,true,true,true,true,true,true,true,true,true,false,
          true,false,true,false,true,false,true,false,true,false,true,true,true,true,true,true,true,true,false,true,false,true,false,true,false,true,false,
          true,false,true,true,true,true,true};
          parameter Modelica.Units.SI.Position cm[n] = {0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,0.01172,0,0.01465,0,0.01758,0,0.02051,0,0.0219875,
          0.0014525,0.0234525,0.0029175,0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,0.010255,0,0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,
          0.002185,0.024185,0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,0.01172,0,0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,0.0029175,
          0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,0.010255,0,0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,0.002185,0.024185};

    parameter Modelica.Units.SI.Length lg[n] = {0.0073,0.00437,0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.011695,0.002905,0.008765,
    0.005835,0.005835,0.008765,0.002905,0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.01316,0.00144,0.01023,0.00437,0.0073,0.0073,0.00437,0.01023,0.00144,
    0.01316,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.011695,0.002905,0.008765,0.005835,0.005835,0.008765,0.002905,0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,
    0,0.0146,0,0.01316,0.00144,0.01023,0.00437,0.0073};

    parameter Modelica.Units.SI.Angle a = 0;
    parameter Modelica.Units.SI.Angle b = 0;
    parameter Modelica.Units.SI.Angle g = 0;

    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b[n] annotation (
      Placement(visible = true, transformation(origin = {-100, -3}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-100, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a[n] annotation (
      Placement(visible = true, transformation(origin = {100, -9}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {100, -8}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    AGILEX_BUNKER_PRO.element_grouser link[n](
      startX=startX,
      startY=startY,
      startZ=startZ,
      thZ=thZ,
      there_is_grouser=there_is_grouser,
      cm=cm,
      lg=lg,
      dh=dh,
      each l=l,
      each a=a,
      each b=b,
      each g=g,
      each start=start,
      each anim=anim) annotation (Placement(visible=true, transformation(
          origin={-4,24},
          extent={{-10,10},{10,-10}},
          rotation=-90)));

  equation
    for n in 1:n-1 loop
      connect(link[n].frame_b, link[n + 1].frame_a);
    end for;
    connect(link[n].frame_b, link[1].frame_a);
    for n in 1:n loop
      connect(link[n].frame_b1, frame_b[n]);
      connect(link[n].frame_a1, frame_a[n]);
    end for;
    annotation (
      Icon(coordinateSystem(grid = {2, 0})),
      Diagram(coordinateSystem(grid = {2, 3})),
      Icon(coordinateSystem(grid = {2, 0})),
      Diagram(coordinateSystem(grid = {2, 3})));
  end track_with_grouser;

  model trackterrainfinal_grouser
   parameter Real start[3] = {0, 0, 0};
   input Modelica.Units.SI.Length dh[80];
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder(diameter = 0.16,
      length=0.105,
      lengthDirection={0,0,1},                                                                                               r = {0, 0, 0.0525}, r_0(each fixed = true, start = start), useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-118, 278}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute(animation = false, useAxisFlange = true) annotation (
      Placement(visible = true, transformation(origin = {-186, 270}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder4(
      r_shape={0,0,-0.0375*2},
      length=0.15,                                                diameter = 0.088,
      r={0,0,0},                                                                                                        useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-124, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute1(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-167, 55}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder2(
      r_shape={0,0,-0.0375*2},
      length=0.15,                                                diameter = 0.088,
      r={0,0,0},                                                                                                        useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin={-126,-8},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute2(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-155, 13}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder6(
      r_shape={0,0,-0.0375*2},
      length=0.15,                                                diameter = 0.088,
      r={0,0,0},                                                                                                        useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-132, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute3(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-173, -67}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder7(
      r_shape={0,0,-0.0375*2},
      length=0.15,                                                diameter = 0.130,
      r={0,0,0},                                                                                                        useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-124, 114}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute4(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-169, 127}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation(animation = true, r = {-0.153, 0.040, 0}) annotation (
      Placement(visible = true, transformation(origin = {-198, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation1(r = {-0.350, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-188, 22}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation2(r = {-0.25, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-190, -6}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute5(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-261, -125}, extent = {{-13, 13}, {13, -13}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation3(animation = false, r = {0.250, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-300, -126}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation4(r = {-0.192, -0.192, 0}) annotation (
      Placement(visible = true, transformation(origin = {-218, -124}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute6(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-181, -121}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder10(
      r_shape={0,0,-0.0375*2},
      length=0.15,                                                 diameter = 0.105,
      r={0,0,0},                                                                                                         useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-128, -148}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute7(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-185, -203}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation5(animation = false, r = {0.281, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-310, -222}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder12(
      r_shape={0,0,-0.0375*2},
      length=0.15,                                                 diameter = 0.105,
      r={0,0,0},                                                                                                          useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-134, -182}, extent = {{-12, -12}, {12, 12}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute8(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-276, -224}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation7(r = {-0.087, -0.029, 0}) annotation (
      Placement(visible = true, transformation(origin = {-286, 212}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder13(
      r_shape={0,0,-0.0375*2},
      length=0.15,                                                 diameter = 0.105,
      r={0,0,0},                                                                                                         useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-122, 194}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute9(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-175, 217}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation9(r = {-0.071, -0.120, 0}) annotation (
      Placement(visible = true, transformation(origin = {-208, 216}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute10(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-264, 240}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation6(r = {-0.192, -0.192, 0}) annotation (
      Placement(visible = true, transformation(origin = {-226, -202}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation8(animation = false, r = {-0.134, -0.134, 0}) annotation (
      Placement(visible = true, transformation(origin = {-288, -84}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation10(animation = false, r = {0.064, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-208, -48}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel springDamperParallel(
      c=100000,
      d=1000,                                                                                             fixedRotationAtFrame_a = false, fixedRotationAtFrame_b = false,
      s_unstretched=0.15)                                                                                                                                                                       annotation (
      Placement(visible = true, transformation(origin = {-264, -58}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation11(animation = false, r = {-0.134, -0.134, 0}) annotation (
      Placement(visible = true, transformation(origin = {-388, -182}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation12(animation = false, r = {0.098, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-376, -10}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel springDamperParallel1(
      c=100000,
      d=1000,                                                                                              fixedRotationAtFrame_a = false, fixedRotationAtFrame_b = false,
      s_unstretched=0.14)                                                                                                                                                                        annotation (
      Placement(visible = true, transformation(origin = {-401, -83}, extent = {{-11, -11}, {11, 11}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel springDamperParallel2(
      c=100000,
      d=1000,                                                                                              fixedRotationAtFrame_a = false, fixedRotationAtFrame_b = false,
      s_unstretched=0.17)                                                                                                                                                                        annotation (
      Placement(visible = true, transformation(origin = {-410, 56}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation13(animation = false, r = {-0.044, -0.065, 0}) annotation (
      Placement(visible = true, transformation(origin = {-250, 142}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation14(animation = false, r = {0.278, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-350, 30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    AGILEX_BUNKER_PRO.contact_sprocket contact_sprocket annotation (Placement(
          visible=true, transformation(
          origin={28,238},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    AGILEX_BUNKER_PRO.track_with_grouser cinghiafinal(start=start, dh=dh)
      annotation (Placement(visible=true, transformation(
          origin={271,67},
          extent={{-67,-67},{67,67}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact(R=0.0525) annotation (Placement(visible=
            true, transformation(
          origin={12,174},
          extent={{-40,-40},{40,40}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact1(R=0.065) annotation (Placement(visible=
            true, transformation(
          origin={6,78},
          extent={{-38,-38},{38,38}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact2(R=0.044) annotation (Placement(visible=
            true, transformation(
          origin={8,-10},
          extent={{-36,-36},{36,36}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact3(R=0.044) annotation (Placement(visible=
            true, transformation(
          origin={-10,-104},
          extent={{-44,-44},{44,44}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact4(R=0.044) annotation (Placement(visible=
            true, transformation(
          origin={-21,-193},
          extent={{-39,-39},{39,39}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact5(R=0.0525) annotation (Placement(visible=
            true, transformation(
          origin={-18,-288},
          extent={{-44,-44},{44,44}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact6(R=0.0525) annotation (Placement(visible=
            true, transformation(
          origin={-129,-347},
          extent={{-49,-49},{49,49}},
          rotation=0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation15(r={-0.105,
          -0.105,-0.0525})                                                                              annotation (
      Placement(visible = true, transformation(origin = {-562, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation16(animation = false, r={0.458,
          0.065,0.0525})                                                                                                annotation (
      Placement(visible = true, transformation(origin = {-528, -44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a annotation (
      Placement(visible = true, transformation(origin = {-664, -36}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-100, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b[80] annotation (
      Placement(visible = true, transformation(origin = {494, 64}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {98, 0}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Components.Damper damper(d=7)
      annotation (Placement(transformation(extent={{-198,306},{-178,326}})));
    Modelica.Mechanics.MultiBody.Joints.Prismatic prismatic(useAxisFlange=true)
                       annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-204,88})));
    Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_a annotation (
        Placement(transformation(extent={{-6,90},{14,110}}), iconTransformation(
            extent={{-6,90},{14,110}})));
    Modelica.Mechanics.Translational.Sources.Position position(exact=true)
      annotation (Placement(transformation(extent={{-252,98},{-232,118}})));
    Modelica.Blocks.Sources.RealExpression realExpression
      annotation (Placement(transformation(extent={{-278,72},{-258,92}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation17(r={
          0.01081078,-0.05765764,0})
      annotation (Placement(transformation(extent={{-254,174},{-234,194}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute11(animation=false)
      annotation (Placement(transformation(extent={{-226,178},{-206,198}})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder16(
      r={0,0,0},
      r_shape={0,0,-0.0375*2},
      length=0.15,
      diameter=0.0756757,
      useQuaternions=false)
      annotation (Placement(transformation(extent={{-188,156},{-168,176}})));
    AGILEX_BUNKER_PRO.contact contact7(R=0.0756757*0.5)
      annotation (Placement(transformation(extent={{-100,162},{-80,182}})));
  equation
    connect(revolute.frame_b, bodyCylinder.frame_a) annotation (
      Line(points = {{-176, 270}, {-140, 270}, {-140, 278}, {-128, 278}}));
    connect(bodyCylinder4.frame_a, revolute1.frame_b) annotation (
      Line(points = {{-134, 60}, {-134, 52.75}, {-138, 52.75}, {-138, 55.5}, {-156, 55.5}, {-156, 55}}, color = {95, 95, 95}));
    connect(bodyCylinder2.frame_a, revolute2.frame_b) annotation (
      Line(points={{-136,-8},{-136,-7.5},{-146,-7.5},{-146,13}},          color = {95, 95, 95}));
    connect(bodyCylinder6.frame_a, revolute3.frame_b) annotation (
      Line(points = {{-142, -70}, {-142, -75.5}, {-162, -75.5}, {-162, -67}}, color = {95, 95, 95}));
    connect(bodyCylinder7.frame_a, revolute4.frame_b) annotation (
      Line(points = {{-134, 114}, {-134, 110.375}, {-158, 110.375}, {-158, 127}}, color = {95, 95, 95}));
    connect(fixedTranslation.frame_b, revolute1.frame_a) annotation (
      Line(points = {{-188, 46}, {-188, 55}, {-178, 55}}));
    connect(fixedTranslation.frame_b, fixedTranslation1.frame_a) annotation (
      Line(points = {{-188, 46}, {-188, 32}}));
    connect(fixedTranslation1.frame_b, revolute2.frame_a) annotation (
      Line(points = {{-188, 12}, {-188, 13}, {-164, 13}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_b, revolute3.frame_a) annotation (
      Line(points = {{-190, -16}, {-191.25, -16}, {-191.25, -32}, {-192.5, -32}, {-192.5, -40}, {-184, -40}, {-184, -67}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_a, fixedTranslation1.frame_b) annotation (
      Line(points = {{-190, 4}, {-189, 4}, {-189, 12}, {-188, 12}}));
    connect(fixedTranslation3.frame_b, revolute5.frame_a) annotation (
      Line(points = {{-290, -126}, {-303, -126}, {-303, -125}, {-274, -125}}, color = {95, 95, 95}));
    connect(revolute5.frame_b, fixedTranslation4.frame_a) annotation (
      Line(points = {{-248, -125}, {-248, -137.5}, {-228, -137.5}, {-228, -124}}));
    connect(bodyCylinder10.frame_a, revolute6.frame_b) annotation (
      Line(points = {{-138, -148}, {-138, -151.813}, {-146, -151.813}, {-146, -149.625}, {-170, -149.625}, {-170, -121}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_b, fixedTranslation3.frame_a) annotation (
      Line(points = {{-190, -16}, {-190, -20}, {-310, -20}, {-310, -126}}, color = {95, 95, 95}));
    connect(bodyCylinder12.frame_a, revolute7.frame_b) annotation (
      Line(points = {{-146, -182}, {-146, -203.5}, {-174, -203.5}, {-174, -203}}, color = {95, 95, 95}));
    connect(fixedTranslation5.frame_b, revolute8.frame_a) annotation (
      Line(points = {{-300, -222}, {-300, -224}, {-286, -224}}, color = {95, 95, 95}));
    connect(fixedTranslation1.frame_b, fixedTranslation5.frame_a) annotation (
      Line(points = {{-188, 12}, {-320, 12}, {-320, -222}}, color = {95, 95, 95}));
    connect(bodyCylinder13.frame_a, revolute9.frame_b) annotation (
      Line(points = {{-132, 194}, {-132, 185.188}, {-144, 185.188}, {-144, 184.375}, {-164, 184.375}, {-164, 217}}, color = {95, 95, 95}));
    connect(revolute9.frame_a, fixedTranslation9.frame_b) annotation (
      Line(points={{-186,217},{-198,217},{-198,216}},        color = {95, 95, 95}));
    connect(revolute10.frame_b, fixedTranslation9.frame_a) annotation (
      Line(points = {{-254, 240}, {-229, 240}, {-229, 216}, {-218, 216}}));
    connect(fixedTranslation7.frame_b, revolute10.frame_a) annotation (
      Line(points = {{-276, 212}, {-274, 212}, {-274, 240}}, color = {95, 95, 95}));
    connect(revolute8.frame_b, fixedTranslation6.frame_a) annotation (
      Line(points = {{-266, -224}, {-240.5, -224}, {-240.5, -202}, {-236, -202}}));
    connect(revolute7.frame_a, fixedTranslation6.frame_b) annotation (
      Line(points = {{-196, -203}, {-216, -203}, {-216, -202}}, color = {95, 95, 95}));
    connect(revolute6.frame_a, fixedTranslation4.frame_b) annotation (
      Line(points = {{-192, -121}, {-208, -121}, {-208, -124}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_b, fixedTranslation10.frame_a) annotation (
      Line(points = {{-190, -16}, {-190, -32}, {-198, -32}, {-198, -48}}, color = {95, 95, 95}));
    connect(fixedTranslation8.frame_a, revolute5.frame_b) annotation (
      Line(points={{-298,-84},{-298,-114},{-248,-114},{-248,-125}}));
    connect(springDamperParallel.frame_b, fixedTranslation8.frame_b) annotation (
      Line(points = {{-264, -68}, {-252, -68}, {-252, -84}, {-278, -84}}, color = {95, 95, 95}));
    connect(springDamperParallel.frame_a, fixedTranslation10.frame_b) annotation (
      Line(points = {{-264, -48}, {-218, -48}}, color = {95, 95, 95}));
    connect(fixedTranslation12.frame_a, fixedTranslation1.frame_b) annotation (
      Line(points = {{-366, -10}, {-220, -10}, {-220, 12}, {-188, 12}}));
    connect(fixedTranslation12.frame_b, springDamperParallel1.frame_a) annotation (
      Line(points = {{-386, -10}, {-401, -10}, {-401, -72}}, color = {95, 95, 95}));
    connect(fixedTranslation11.frame_a, revolute8.frame_b) annotation (
      Line(points = {{-378, -182}, {-266, -182}, {-266, -224}}, color = {95, 95, 95}));
    connect(fixedTranslation11.frame_b, springDamperParallel1.frame_b) annotation (
      Line(points={{-398,-182},{-401,-182},{-401,-94}},        color = {95, 95, 95}));
    connect(fixedTranslation14.frame_a, fixedTranslation1.frame_b) annotation (
      Line(points = {{-340, 30}, {-188, 30}, {-188, 12}}, color = {95, 95, 95}));
    connect(fixedTranslation14.frame_b, springDamperParallel2.frame_a) annotation (
      Line(points = {{-360, 30}, {-410, 30}, {-410, 46}}, color = {95, 95, 95}));
    connect(revolute10.frame_b, fixedTranslation13.frame_a) annotation (
      Line(points = {{-254, 240}, {-224, 240}, {-224, 152}, {-250, 152}}));
    connect(fixedTranslation13.frame_b, springDamperParallel2.frame_b) annotation (
      Line(points = {{-250, 132}, {-410, 132}, {-410, 66}}));
    connect(bodyCylinder.frame_b, contact_sprocket.frame_b) annotation (
      Line(points={{-108,278},{16.2,278},{16.2,238.4}}));
    connect(contact.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{52,174.8},{204,174.8},{204,65.66}}, color = {95, 95, 95}, thickness = 0.5));
    connect(contact1.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{44,78.76},{204,78.76},{204,65.66}},
                                                      thickness = 0.5));
    connect(contact2.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{44,-9.28},{204,-9.28},{204,65.66}}, thickness = 0.5));
    connect(contact3.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{34,-103.12},{204,-103.12},{204,65.66}},
                                                          color = {95, 95, 95}, thickness = 0.5));
    connect(contact4.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{18,-192.22},{204,-192.22},{204,65.66}},
                                                          thickness = 0.5));
    connect(contact5.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{26,-287.12},{204,-287.12},{204,65.66}},
                                                          color = {95, 95, 95}, thickness = 0.5));
    connect(contact6.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{-80,-346.02},{204,-346.02},{204,65.66}},
                                                           thickness = 0.5));
    connect(contact_sprocket.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{38,238.2},{124,238.2},{124,65.66},{204,65.66}},thickness = 0.5));
    connect(fixedTranslation15.frame_a, fixedTranslation2.frame_b) annotation (
      Line(points = {{-562, 42}, {-586, 42}, {-586, -16}, {-190, -16}}));
    connect(fixedTranslation15.frame_b, revolute.frame_a) annotation (
      Line(points = {{-562, 62}, {-448, 62}, {-448, 270}, {-196, 270}}, color = {95, 95, 95}));
    connect(fixedTranslation.frame_a, fixedTranslation7.frame_a) annotation (
      Line(points = {{-208, 46}, {-320, 46}, {-320, 212}, {-296, 212}}, color = {95, 95, 95}));
    connect(fixedTranslation16.frame_b, fixedTranslation.frame_a) annotation (
      Line(points = {{-518, -44}, {-208, -44}, {-208, 46}}));
    connect(frame_a, fixedTranslation16.frame_a) annotation (
      Line(points = {{-664, -36}, {-538, -36}, {-538, -44}}));
    connect(cinghiafinal.frame_a, frame_b) annotation (
      Line(points={{338,61.64},{494,61.64},{494,64}},  thickness = 0.5));
    connect(damper.flange_b, revolute.axis) annotation (Line(points={{-178,316},
            {-162,316},{-162,280},{-186,280}}, color={0,0,0}));
    connect(damper.flange_a, revolute.support) annotation (Line(points={{-198,
            316},{-216,316},{-216,300},{-222,300},{-222,280},{-192,280}}, color=
           {0,0,0}));
    connect(prismatic.frame_b, revolute4.frame_a) annotation (Line(
        points={{-204,98},{-204,130},{-180,130},{-180,127}},
        color={95,95,95},
        thickness=0.5));
    connect(prismatic.frame_a, fixedTranslation.frame_a) annotation (Line(
        points={{-204,78},{-204,58},{-208,58},{-208,46}},
        color={95,95,95},
        thickness=0.5));
    connect(flange_a, revolute.axis) annotation (Line(points={{4,100},{-280,100},
            {-280,330},{-186,330},{-186,280}}, color={0,0,0}));
    connect(position.flange, prismatic.axis) annotation (Line(points={{-232,108},{
            -226,108},{-226,94},{-210,94},{-210,96}}, color={0,127,0}));
    connect(realExpression.y, position.s_ref) annotation (Line(points={{-257,82},{
            -256,82},{-256,108},{-254,108}}, color={0,0,127}));
    connect(revolute10.frame_b, fixedTranslation17.frame_a) annotation (Line(
        points={{-254,240},{-254,184}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute11.frame_a, fixedTranslation17.frame_b) annotation (Line(
        points={{-226,188},{-228,188},{-228,184},{-234,184}},
        color={95,95,95},
        thickness=0.5));
    connect(contact7.frame_a, cinghiafinal.frame_b) annotation (Line(
        points={{-80,172.2},{-80,58},{-36,58},{-36,36},{118,36},{118,65.66},{
            204,65.66}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute11.frame_b, bodyCylinder16.frame_a) annotation (Line(
        points={{-206,188},{-206,166},{-188,166},{-188,166}},
        color={95,95,95},
        thickness=0.5));
    connect(contact.frame_b, bodyCylinder13.frame_b) annotation (Line(
        points={{-35.2,175.6},{-35.2,194},{-112,194}},
        color={95,95,95},
        thickness=0.5));
    connect(contact7.frame_b, bodyCylinder16.frame_b) annotation (Line(
        points={{-101.8,172.4},{-148,172.4},{-148,166},{-168,166}},
        color={95,95,95},
        thickness=0.5));
    connect(contact1.frame_b, bodyCylinder7.frame_b) annotation (Line(
        points={{-38.84,79.52},{-72,79.52},{-72,114},{-114,114}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder4.frame_b, contact2.frame_b) annotation (Line(
        points={{-114,60},{-94,60},{-94,28},{-34.48,28},{-34.48,-8.56}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder2.frame_b, contact3.frame_b) annotation (Line(
        points={{-116,-8},{-90,-8},{-90,-26},{-61.92,-26},{-61.92,-102.24}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder6.frame_b, contact4.frame_b) annotation (Line(
        points={{-122,-70},{-102,-70},{-102,-74},{-74,-74},{-74,-191.44},{
            -67.02,-191.44}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder10.frame_b, contact5.frame_b) annotation (Line(
        points={{-118,-148},{-118,-166},{-100,-166},{-100,-226},{-69.92,-226},{
            -69.92,-286.24}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder12.frame_b, contact6.frame_b) annotation (Line(
        points={{-122,-182},{-122,-286},{-208,-286},{-208,-345.04},{-186.82,
            -345.04}},
        color={95,95,95},
        thickness=0.5));
    annotation (
      Icon(coordinateSystem(grid = {2, 0})),
      experiment(StartTime = 0, StopTime = 10, Tolerance = 0.0001, Interval = 0.02));
  end trackterrainfinal_grouser;

  model trackterrainfinal1_grouser
   parameter Real start[3] = {0, 0, 0};
   input Modelica.Units.SI.Length dh[80];
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder(diameter = 0.16,
      length=0.105,
      lengthDirection={0,0,1},                                                                                               r = {0, 0, 0.0525}, useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-120, 276}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute(animation = false, useAxisFlange = true) annotation (
      Placement(visible = true, transformation(origin = {-186, 270}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder4(
      r_shape={0,0,-0.0375*2},
      length=0.15,                                                diameter = 0.088,
      r={0,0,0},                                                                                                        useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-124, 60}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute1(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-167, 55}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder2(
      r_shape={0,0,-0.0375*2},
      length=0.15,                                                diameter = 0.088,
      r={0,0,0},                                                                                                        useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-125, -3}, extent = {{-9, -9}, {9, 9}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute2(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-157, 11}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder6(
      r_shape={0,0,-0.0375*2},
      length=0.15,                                                diameter = 0.088,
      r={0,0,0},                                                                                                        useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-132, -70}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute3(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-173, -67}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder7(
      r_shape={0,0,-0.0375*2},
      length=0.15,                                                diameter = 0.130,
      r={0,0,0},                                                                                                        useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-124, 114}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute4(animation = false) annotation (
      Placement(visible = true, transformation(origin={-167,127},    extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation(animation = true, r={-0.153,
          0.040,0})                                                                                                annotation (
      Placement(visible = true, transformation(origin={-200,46},    extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation1(r = {-0.350, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-188, 22}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation2(r = {-0.25, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-190, -6}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute5(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-261, -125}, extent = {{-13, 13}, {13, -13}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation3(animation = false, r = {0.250, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-300, -126}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation4(r = {-0.192, -0.192, 0}) annotation (
      Placement(visible = true, transformation(origin = {-218, -124}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute6(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-181, -121}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder10(
      r_shape={0,0,-0.0375*2},
      length=0.15,                                                 diameter = 0.105,
      r={0,0,0},                                                                                                         useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-128, -148}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder11(
      r_shape={0,0,-0.0375*2},
      length=0.15,                                                 diameter = 0.105,
      r={0,0,0},                                                                                                         useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-134, -222}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute7(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-185, -203}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation5(animation = false, r = {0.281, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-310, -222}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute8(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-276, -224}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation7(r = {-0.087, -0.029, 0}) annotation (
      Placement(visible = true, transformation(origin = {-288, 212}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder13(
      r_shape={0,0,-0.0375*2},
      length=0.15,                                                 diameter = 0.105,
      r={0,0,0},                                                                                                         useQuaternions = false) annotation (
      Placement(visible = true, transformation(origin = {-122, 194}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute9(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-175, 217}, extent = {{-11, -11}, {11, 11}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation9(r = {-0.071, -0.120, 0}) annotation (
      Placement(visible = true, transformation(origin = {-208, 216}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute10(animation = false) annotation (
      Placement(visible = true, transformation(origin = {-246, 218}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation6(r = {-0.192, -0.192, 0}) annotation (
      Placement(visible = true, transformation(origin = {-226, -202}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation8(animation = false, r = {-0.134, -0.134, 0}) annotation (
      Placement(visible = true, transformation(origin = {-288, -84}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation10(animation = false, r = {0.064, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-208, -48}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel springDamperParallel(
      c=100000,
      d=1000,                                                                                             fixedRotationAtFrame_a = false, fixedRotationAtFrame_b = false,
      s_unstretched=0.15)                                                                                                                                                                       annotation (
      Placement(visible = true, transformation(origin = {-264, -58}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation11(animation = false, r = {-0.134, -0.134, 0}) annotation (
      Placement(visible = true, transformation(origin = {-388, -182}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation12(animation = false, r = {0.098, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-376, -10}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel springDamperParallel1(
      c=100000,
      d=1000,                                                                                              fixedRotationAtFrame_a = false, fixedRotationAtFrame_b = false,
      s_unstretched=0.14)                                                                                                                                                                        annotation (
      Placement(visible = true, transformation(origin = {-401, -83}, extent = {{-11, -11}, {11, 11}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel springDamperParallel2(
      c=100000,
      d=1000,                                                                                              fixedRotationAtFrame_a = false, fixedRotationAtFrame_b = false,
      s_unstretched=0.17)                                                                                                                                                                        annotation (
      Placement(visible = true, transformation(origin = {-410, 56}, extent = {{10, -10}, {-10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation13(animation = false, r = {-0.044, -0.065, 0}) annotation (
      Placement(visible = true, transformation(origin = {-250, 142}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation14(animation = false, r = {0.278, 0, 0}) annotation (
      Placement(visible = true, transformation(origin = {-350, 30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    AGILEX_BUNKER_PRO.contact_sprocket contact_sprocket annotation (Placement(
          visible=true, transformation(
          origin={28,238},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    AGILEX_BUNKER_PRO.track_with_grouser cinghiafinal(start=start, dh=dh)
      annotation (Placement(visible=true, transformation(
          origin={271,67},
          extent={{-67,-67},{67,67}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact(R=0.0525) annotation (Placement(visible=
            true, transformation(
          origin={12,174},
          extent={{-40,-40},{40,40}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact1(R=0.065) annotation (Placement(visible=
            true, transformation(
          origin={6,78},
          extent={{-38,-38},{38,38}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact2(R=0.044) annotation (Placement(visible=
            true, transformation(
          origin={8,-10},
          extent={{-36,-36},{36,36}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact3(R=0.044) annotation (Placement(visible=
            true, transformation(
          origin={-10,-104},
          extent={{-44,-44},{44,44}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact4(R=0.044) annotation (Placement(visible=
            true, transformation(
          origin={-21,-193},
          extent={{-39,-39},{39,39}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact5(R=0.0525) annotation (Placement(visible=
            true, transformation(
          origin={-18,-288},
          extent={{-44,-44},{44,44}},
          rotation=0)));
    AGILEX_BUNKER_PRO.contact contact6(R=0.0525) annotation (Placement(visible=
            true, transformation(
          origin={-129,-347},
          extent={{-49,-49},{49,49}},
          rotation=0)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation15(r={-0.105,
          -0.105,-0.0525})                                                                              annotation (
      Placement(visible = true, transformation(origin = {-562, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation16(animation = false, r={0.458,
          0.065,-0.0525})                                                                                                annotation (
      Placement(visible = true, transformation(origin = {-530, -46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a annotation (
      Placement(visible = true, transformation(origin = {-664, -36}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-100, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b[80] annotation (
      Placement(visible = true, transformation(origin = {494, 64}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {98, 0}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Mechanics.Rotational.Components.Damper damper(d=7)
      annotation (Placement(transformation(extent={{-200,324},{-180,344}})));
    Modelica.Mechanics.MultiBody.Joints.Prismatic prismatic(useAxisFlange=true)
                       annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-202,96})));
    Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_a annotation (
        Placement(transformation(extent={{-6,86},{14,106}}), iconTransformation(
            extent={{-6,86},{14,106}})));
    Modelica.Mechanics.Translational.Sources.Position position(exact=true)
      annotation (Placement(transformation(extent={{-248,84},{-228,104}})));
    Modelica.Blocks.Sources.RealExpression realExpression
      annotation (Placement(transformation(extent={{-278,72},{-258,92}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation17(r={
          0.01081078,-0.05765764,0})
      annotation (Placement(transformation(extent={{-240,168},{-220,188}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute11(animation=false)
      annotation (Placement(transformation(extent={{-212,168},{-192,188}})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder16(
      r={0,0,0},
      r_shape={0,0,-0.0375*2},
      length=0.15,
      diameter=0.0756757,
      useQuaternions=false)
      annotation (Placement(transformation(extent={{-180,156},{-160,176}})));
    AGILEX_BUNKER_PRO.contact contact7(R=0.0756757*0.5)
      annotation (Placement(transformation(extent={{-102,146},{-82,166}})));
  equation
    connect(revolute.frame_b, bodyCylinder.frame_a) annotation (
      Line(points = {{-176, 270}, {-140, 270}, {-140, 276}, {-130, 276}}));
    connect(bodyCylinder4.frame_a, revolute1.frame_b) annotation (
      Line(points = {{-134, 60}, {-134, 52.75}, {-138, 52.75}, {-138, 55.5}, {-156, 55.5}, {-156, 55}}, color = {95, 95, 95}));
    connect(bodyCylinder2.frame_a, revolute2.frame_b) annotation (
      Line(points = {{-134, -3}, {-134, -7.5}, {-146, -7.5}, {-146, 11}}, color = {95, 95, 95}));
    connect(bodyCylinder6.frame_a, revolute3.frame_b) annotation (
      Line(points = {{-142, -70}, {-142, -75.5}, {-162, -75.5}, {-162, -67}}, color = {95, 95, 95}));
    connect(bodyCylinder7.frame_a, revolute4.frame_b) annotation (
      Line(points={{-134,114},{-134,110.375},{-156,110.375},{-156,127}},          color = {95, 95, 95}));
    connect(fixedTranslation.frame_b, revolute1.frame_a) annotation (
      Line(points={{-190,46},{-190,55},{-178,55}}));
    connect(fixedTranslation.frame_b, fixedTranslation1.frame_a) annotation (
      Line(points={{-190,46},{-190,40},{-188,40},{-188,32}}));
    connect(fixedTranslation1.frame_b, revolute2.frame_a) annotation (
      Line(points = {{-188, 12}, {-188, 11}, {-168, 11}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_b, revolute3.frame_a) annotation (
      Line(points = {{-190, -16}, {-191.25, -16}, {-191.25, -32}, {-192.5, -32}, {-192.5, -40}, {-184, -40}, {-184, -67}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_a, fixedTranslation1.frame_b) annotation (
      Line(points = {{-190, 4}, {-189, 4}, {-189, 12}, {-188, 12}}));
    connect(fixedTranslation3.frame_b, revolute5.frame_a) annotation (
      Line(points = {{-290, -126}, {-303, -126}, {-303, -125}, {-274, -125}}, color = {95, 95, 95}));
    connect(revolute5.frame_b, fixedTranslation4.frame_a) annotation (
      Line(points = {{-248, -125}, {-248, -137.5}, {-228, -137.5}, {-228, -124}}));
    connect(bodyCylinder10.frame_a, revolute6.frame_b) annotation (
      Line(points = {{-138, -148}, {-138, -151.813}, {-146, -151.813}, {-146, -149.625}, {-170, -149.625}, {-170, -121}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_b, fixedTranslation3.frame_a) annotation (
      Line(points = {{-190, -16}, {-190, -20}, {-310, -20}, {-310, -126}}, color = {95, 95, 95}));
    connect(bodyCylinder11.frame_a, revolute7.frame_b) annotation (
      Line(points = {{-144, -222}, {-144, -238.813}, {-174, -238.813}, {-174, -203}}, color = {95, 95, 95}));
    connect(fixedTranslation5.frame_b, revolute8.frame_a) annotation (
      Line(points = {{-300, -222}, {-300, -224}, {-286, -224}}, color = {95, 95, 95}));
    connect(fixedTranslation1.frame_b, fixedTranslation5.frame_a) annotation (
      Line(points = {{-188, 12}, {-320, 12}, {-320, -222}}, color = {95, 95, 95}));
    connect(bodyCylinder13.frame_a, revolute9.frame_b) annotation (
      Line(points = {{-132, 194}, {-132, 185.188}, {-144, 185.188}, {-144, 184.375}, {-164, 184.375}, {-164, 217}}, color = {95, 95, 95}));
    connect(revolute9.frame_a, fixedTranslation9.frame_b) annotation (
      Line(points={{-186,217},{-198,217},{-198,216}},        color = {95, 95, 95}));
    connect(revolute10.frame_b, fixedTranslation9.frame_a) annotation (
      Line(points = {{-236, 218}, {-223, 218}, {-223, 216}, {-218, 216}}));
    connect(revolute8.frame_b, fixedTranslation6.frame_a) annotation (
      Line(points = {{-266, -224}, {-240.5, -224}, {-240.5, -202}, {-236, -202}}));
    connect(revolute7.frame_a, fixedTranslation6.frame_b) annotation (
      Line(points = {{-196, -203}, {-216, -203}, {-216, -202}}, color = {95, 95, 95}));
    connect(revolute6.frame_a, fixedTranslation4.frame_b) annotation (
      Line(points = {{-192, -121}, {-208, -121}, {-208, -124}}, color = {95, 95, 95}));
    connect(fixedTranslation2.frame_b, fixedTranslation10.frame_a) annotation (
      Line(points = {{-190, -16}, {-190, -32}, {-198, -32}, {-198, -48}}, color = {95, 95, 95}));
    connect(fixedTranslation8.frame_a, revolute5.frame_b) annotation (
      Line(points={{-298,-84},{-298,-114},{-248,-114},{-248,-125}}));
    connect(springDamperParallel.frame_b, fixedTranslation8.frame_b) annotation (
      Line(points = {{-264, -68}, {-252, -68}, {-252, -84}, {-278, -84}}, color = {95, 95, 95}));
    connect(springDamperParallel.frame_a, fixedTranslation10.frame_b) annotation (
      Line(points = {{-264, -48}, {-218, -48}}, color = {95, 95, 95}));
    connect(fixedTranslation12.frame_a, fixedTranslation1.frame_b) annotation (
      Line(points = {{-366, -10}, {-220, -10}, {-220, 12}, {-188, 12}}));
    connect(fixedTranslation12.frame_b, springDamperParallel1.frame_a) annotation (
      Line(points = {{-386, -10}, {-401, -10}, {-401, -72}}, color = {95, 95, 95}));
    connect(fixedTranslation11.frame_a, revolute8.frame_b) annotation (
      Line(points = {{-378, -182}, {-266, -182}, {-266, -224}}, color = {95, 95, 95}));
    connect(fixedTranslation11.frame_b, springDamperParallel1.frame_b) annotation (
      Line(points = {{-398, -182}, {-401, -182}, {-401, -94}}, color = {95, 95, 95}));
    connect(fixedTranslation14.frame_a, fixedTranslation1.frame_b) annotation (
      Line(points={{-340,30},{-188,30},{-188,12}},        color = {95, 95, 95}));
    connect(fixedTranslation14.frame_b, springDamperParallel2.frame_a) annotation (
      Line(points = {{-360, 30}, {-410, 30}, {-410, 46}}, color = {95, 95, 95}));
    connect(revolute10.frame_b, fixedTranslation13.frame_a) annotation (
      Line(points = {{-236, 218}, {-224, 218}, {-224, 152}, {-250, 152}}));
    connect(fixedTranslation13.frame_b, springDamperParallel2.frame_b) annotation (
      Line(points = {{-250, 132}, {-410, 132}, {-410, 66}}));
    connect(bodyCylinder.frame_b, contact_sprocket.frame_b) annotation (
      Line(points={{-110,276},{16.2,276},{16.2,238.4}}));
    connect(contact.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{52,174.8},{204,174.8},{204,65.66}}, color = {95, 95, 95}, thickness = 0.5));
    connect(contact1.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{44,78.76},{204,78.76},{204,65.66}},
                                                      thickness = 0.5));
    connect(contact2.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{44,-9.28},{204,-9.28},{204,65.66}}, thickness = 0.5));
    connect(contact3.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{34,-103.12},{204,-103.12},{204,65.66}},
                                                          color = {95, 95, 95}, thickness = 0.5));
    connect(contact4.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{18,-192.22},{204,-192.22},{204,65.66}},
                                                          thickness = 0.5));
    connect(contact5.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{26,-287.12},{204,-287.12},{204,65.66}},
                                                          color = {95, 95, 95}, thickness = 0.5));
    connect(contact6.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{-80,-346.02},{204,-346.02},{204,65.66}},
                                                           thickness = 0.5));
    connect(contact_sprocket.frame_a, cinghiafinal.frame_b) annotation (
      Line(points={{38,238.2},{124,238.2},{124,65.66},{204,65.66}},thickness = 0.5));
  connect(fixedTranslation15.frame_a, fixedTranslation2.frame_b) annotation (
      Line(points = {{-562, 42}, {-586, 42}, {-586, -16}, {-190, -16}}));
  connect(fixedTranslation15.frame_b, revolute.frame_a) annotation (
      Line(points = {{-562, 62}, {-448, 62}, {-448, 270}, {-196, 270}}, color = {95, 95, 95}));
    connect(fixedTranslation16.frame_b, fixedTranslation.frame_a) annotation (
      Line(points={{-520,-46},{-210,-46},{-210,46}}));
    connect(frame_a, fixedTranslation16.frame_a) annotation (
      Line(points = {{-664, -36}, {-540, -36}, {-540, -46}}));
    connect(cinghiafinal.frame_a, frame_b) annotation (
      Line(points={{338,61.64},{494,61.64},{494,64}},  thickness = 0.5));
    connect(fixedTranslation7.frame_b, revolute10.frame_a) annotation (
      Line(points = {{-278, 212}, {-256, 212}, {-256, 218}}, color = {95, 95, 95}));
    connect(damper.flange_b, revolute.axis) annotation (Line(points={{-180,334},
            {-180,307},{-186,307},{-186,280}}, color={0,0,0}));
    connect(damper.flange_a, revolute.support) annotation (Line(points={{-200,
            334},{-200,307},{-192,307},{-192,280}}, color={0,0,0}));
    connect(fixedTranslation7.frame_a, fixedTranslation.frame_a) annotation (
        Line(
        points={{-298,212},{-358,212},{-358,46},{-210,46}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute4.frame_a, prismatic.frame_b) annotation (Line(
        points={{-178,127},{-202,127},{-202,106}},
        color={95,95,95},
        thickness=0.5));
    connect(prismatic.frame_a, fixedTranslation.frame_a) annotation (Line(
        points={{-202,86},{-202,46},{-210,46}},
        color={95,95,95},
        thickness=0.5));
    connect(flange_a, revolute.axis) annotation (Line(points={{4,96},{-424,96},
            {-424,350},{-186,350},{-186,280}}, color={0,0,0}));
    connect(position.flange, prismatic.axis) annotation (Line(points={{-228,94},{-230,
            94},{-230,104},{-208,104}}, color={0,127,0}));
    connect(realExpression.y, position.s_ref) annotation (Line(points={{-257,82},{
            -252,82},{-252,90},{-250,90},{-250,94}}, color={0,0,127}));
    connect(revolute10.frame_b, fixedTranslation17.frame_a) annotation (Line(
        points={{-236,218},{-236,198},{-240,198},{-240,178}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation17.frame_b, revolute11.frame_a) annotation (Line(
        points={{-220,178},{-212,178}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute11.frame_b, bodyCylinder16.frame_a) annotation (Line(
        points={{-192,178},{-188,178},{-188,164},{-180,164},{-180,166}},
        color={95,95,95},
        thickness=0.5));
    connect(contact7.frame_a, cinghiafinal.frame_b) annotation (Line(
        points={{-82,156.2},{-82,58},{-36,58},{-36,36},{118,36},{118,65.66},{
            204,65.66}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder13.frame_b, contact.frame_b) annotation (Line(
        points={{-112,194},{-96,194},{-96,180},{-35.2,180},{-35.2,175.6}},
        color={95,95,95},
        thickness=0.5));
    connect(contact7.frame_b, bodyCylinder16.frame_b) annotation (Line(
        points={{-103.8,156.4},{-130,156.4},{-130,166},{-160,166}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder7.frame_b, contact1.frame_b) annotation (Line(
        points={{-114,114},{-38.84,114},{-38.84,79.52}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder4.frame_b, contact2.frame_b) annotation (Line(
        points={{-114,60},{-96,60},{-96,32},{-34.48,32},{-34.48,-8.56}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder2.frame_b, contact3.frame_b) annotation (Line(
        points={{-116,-3},{-61.92,-3},{-61.92,-102.24}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder6.frame_b, contact4.frame_b) annotation (Line(
        points={{-122,-70},{-67.02,-70},{-67.02,-191.44}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder10.frame_b, contact5.frame_b) annotation (Line(
        points={{-118,-148},{-108,-148},{-108,-158},{-94,-158},{-94,-222},{
            -69.92,-222},{-69.92,-286.24}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder11.frame_b, contact6.frame_b) annotation (Line(
        points={{-124,-222},{-114,-222},{-114,-230},{-112,-230},{-112,-282},{
            -202,-282},{-202,-345.04},{-186.82,-345.04}},
        color={95,95,95},
        thickness=0.5));
    annotation (
      Icon(coordinateSystem(grid = {2, 0})),
      experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-4, Interval = 0.02));
  end trackterrainfinal1_grouser;

  model track_bunker_pro
    extends Modelica.Icons.Example;
    Modelica.Mechanics.MultiBody.Parts.Fixed fixed
      annotation (Placement(transformation(extent={{-148,-148},{-128,-128}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute(useAxisFlange=true,
        animation=false)
      annotation (Placement(transformation(extent={{-64,-90},{-44,-70}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation(r={-0.1749185,
          -0.0689022,0}) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-50,-134})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute2(useAxisFlange=false,
        animation=false) annotation (Placement(transformation(
          extent={{-10,-11},{10,11}},
          rotation=270,
          origin={-62,-169})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder1(
      r={0,0,0.105*0.5},
      r_shape={0,0,-0.105*0.5},
      diameter=0.16) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-44,-238})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation1(r={
          0.216,-0.042,0}) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-74,-52})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation2(r={
          0.224,0.01961,0}) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-74,-18})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation3(r={-0.033266,
          0.156504,0})
      annotation (Placement(transformation(extent={{-94,-54},{-114,-34}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation4(r={
          0.206,0.029,0}) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-96,66})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation5(r={
          0.104,-0.094,0}) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-66,66})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute1(animation=false)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-94,90})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute3(animation=false)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-64,90})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder(
      r={0,0,0.15*0.5},
      r_shape={0,0,-0.15*0.5},
      diameter=0.135) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-98,122})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder2(
      r={0,0,0.15*0.5},
      r_shape={0,0,-0.15*0.5},
      diameter=0.125) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-56,122})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute4(animation=false)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={-136,-44})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder3(
      r={0,0,0.15*0.5},
      r_shape={0,0,-0.15*0.5},
      diameter=0.08)
      annotation (Placement(transformation(extent={{-156,-54},{-176,-34}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation6(r={
          0.17067,0.088845,0})
      annotation (Placement(transformation(extent={{-102,-28},{-122,-8}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation7(r={
          0.065,-0.076,0})
      annotation (Placement(transformation(extent={{-62,-76},{-42,-56}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute5(useAxisFlange=true,
        animation=false)
      annotation (Placement(transformation(extent={{-48,-50},{-28,-30}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute6(useAxisFlange=true,
        animation=false)
      annotation (Placement(transformation(extent={{-48,-22},{-28,-2}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation8(r={
          0.076368,0.076368,0})
      annotation (Placement(transformation(extent={{-10,-40},{10,-20}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation9(r={
          0.115,-0.090,0})
      annotation (Placement(transformation(extent={{-14,-70},{6,-50}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation10(r={-0.076368,
          0.076368,0})
      annotation (Placement(transformation(extent={{-4,16},{16,36}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation11(r={-0.115,
          -0.090,0})
      annotation (Placement(transformation(extent={{-8,-14},{12,6}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute7(useAxisFlange=true,
        animation=false)
      annotation (Placement(transformation(extent={{24,-14},{44,6}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute8(useAxisFlange=true,
        animation=false)
      annotation (Placement(transformation(extent={{18,-70},{38,-50}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation12(r={
          0.051357,-0.016687,0})
      annotation (Placement(transformation(extent={{72,4},{92,24}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation13(r={-0.051357,
          -0.016687,0})
      annotation (Placement(transformation(extent={{72,-20},{92,0}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation14(r={
          0.051357,-0.016687,0})
      annotation (Placement(transformation(extent={{60,-46},{80,-26}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation15(r={-0.051357,
          -0.016687,0})
      annotation (Placement(transformation(extent={{62,-98},{82,-78}})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder4(
      r={0,0,0.023*0.5},
      r_shape={0,0,-0.023*0.5},
      diameter=0.07)
      annotation (Placement(transformation(extent={{108,4},{128,24}})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder5(
      r={0,0,0.023*0.5},
      r_shape={0,0,-0.023*0.5},
      diameter=0.07)
      annotation (Placement(transformation(extent={{108,-18},{128,2}})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder6(
      r={0,0,0.023*0.5},
      r_shape={0,0,-0.023*0.5},
      diameter=0.07)
      annotation (Placement(transformation(extent={{98,-44},{118,-24}})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder7(
      r={0,0,0.023*0.5},
      r_shape={0,0,-0.023*0.5},
      diameter=0.07)
      annotation (Placement(transformation(extent={{106,-82},{126,-62}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation16(r={
          0.0949,0.0373831,0})
      annotation (Placement(transformation(extent={{-14,-138},{6,-118}})));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel
      springDamperParallel(
      c=1000000,
      s_unstretched=0.12,
      d=10000) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={4,-102})));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel
      springDamperParallel1(
      c=1000000,
      s_unstretched=0.15,
      d=10000) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={48,60})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation17(r={
          0.206,0.029,0} + {-0.12,0.04,0}) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-120,26})));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel
      springDamperParallel2(
      c=1000000,
      s_unstretched=0.14,
      d=10000) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-156,6})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute9(useAxisFlange=true,
        animation=false) annotation (Placement(transformation(
          extent={{-7,-7},{7,7}},
          rotation=90,
          origin={-87,3})));
    Modelica.Mechanics.MultiBody.Joints.Prismatic prismatic(useAxisFlange=true)
                       annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-86,30})));
    Modelica.Mechanics.Translational.Sources.Position position(exact=true)
      annotation (Placement(transformation(extent={{-148,56},{-128,76}})));
    Modelica.Blocks.Sources.RealExpression realExpression
      annotation (Placement(transformation(extent={{-192,56},{-172,76}})));
    cinghia cinghia1(
      sX=1.130312455 - 0.2,
      sZ=0,
      l=0.0586*38)
      annotation (Placement(transformation(extent={{268,-48},{320,4}})));
    contact contact1(R=0.07*0.5)
      annotation (Placement(transformation(extent={{156,20},{176,40}})));
    contact contact2(R=0.07*0.5)
      annotation (Placement(transformation(extent={{156,-18},{176,2}})));
    contact contact3(R=0.07*0.5)
      annotation (Placement(transformation(extent={{152,-50},{172,-30}})));
    contact contact4(R=0.07*0.5)
      annotation (Placement(transformation(extent={{144,-92},{164,-72}})));
    contact contact5(R=0.125*0.5)
      annotation (Placement(transformation(extent={{-12,96},{8,116}})));
    contact contact6(R=0.135*0.5)
      annotation (Placement(transformation(extent={{-44,152},{-24,172}})));
    contact contact7(R=0.08*0.5)
      annotation (Placement(transformation(extent={{-192,-40},{-212,-20}})));
    contact_sprocket contact_sprocket1
      annotation (Placement(transformation(extent={{102,-196},{122,-176}})));
    rotational_stop rotational_stop1
      annotation (Placement(transformation(extent={{30,12},{46,26}})));
    rotational_stop rotational_stop2
      annotation (Placement(transformation(extent={{24,-44},{40,-30}})));
    Modelica.Mechanics.Rotational.Sources.Position position1(exact=true)
      annotation (Placement(transformation(extent={{-156,-86},{-136,-66}})));
    Modelica.Blocks.Sources.RealExpression realExpression1
      annotation (Placement(transformation(extent={{-210,-86},{-190,-66}})));
    Modelica.Mechanics.Rotational.Sources.Position position2(exact=true)
      annotation (Placement(transformation(extent={{-188,34},{-168,54}})));
    Modelica.Blocks.Sources.RealExpression realExpression2
      annotation (Placement(transformation(extent={{-242,34},{-222,54}})));
    Modelica.Mechanics.Rotational.Sources.Position position3(exact=true)
      annotation (Placement(transformation(extent={{-68,14},{-48,34}})));
    Modelica.Blocks.Sources.RealExpression realExpression3
      annotation (Placement(transformation(extent={{-58,32},{-38,52}})));
    Modelica.Mechanics.Rotational.Sources.Position position4(exact=true)
      annotation (Placement(transformation(extent={{-146,-110},{-126,-90}})));
    Modelica.Blocks.Sources.RealExpression realExpression4
      annotation (Placement(transformation(extent={{-220,-116},{-200,-96}})));
    Modelica.Mechanics.Rotational.Sources.Position position5(exact=true)
      annotation (Placement(transformation(extent={{96,-136},{116,-116}})));
    Modelica.Blocks.Sources.RealExpression realExpression5
      annotation (Placement(transformation(extent={{42,-136},{62,-116}})));
    Modelica.Mechanics.Rotational.Sources.Position position6(exact=true)
      annotation (Placement(transformation(extent={{120,44},{140,64}})));
    Modelica.Blocks.Sources.RealExpression realExpression6
      annotation (Placement(transformation(extent={{66,44},{86,64}})));
  equation
    connect(revolute.frame_a, fixed.frame_b) annotation (Line(
        points={{-64,-80},{-74,-80},{-74,-120},{-128,-120},{-128,-138}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute.frame_b, fixedTranslation.frame_a) annotation (Line(
        points={{-44,-80},{-38,-80},{-38,-124},{-50,-124}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute2.frame_a, fixedTranslation.frame_b) annotation (Line(
        points={{-62,-159},{-62,-150},{-50,-150},{-50,-144}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute2.frame_b, bodyCylinder1.frame_a) annotation (Line(
        points={{-62,-179},{-62,-208},{-44,-208},{-44,-228}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation1.frame_a, fixed.frame_b) annotation (Line(
        points={{-74,-62},{-74,-120},{-128,-120},{-128,-138}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation2.frame_a, fixedTranslation1.frame_b) annotation (
        Line(
        points={{-74,-28},{-74,-42}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation3.frame_a, fixedTranslation1.frame_b) annotation (
        Line(
        points={{-94,-44},{-88,-44},{-88,-34},{-74,-34},{-74,-42}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute1.frame_a, fixedTranslation4.frame_b) annotation (Line(
        points={{-94,80},{-94,76},{-96,76}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute3.frame_a, fixedTranslation5.frame_b) annotation (Line(
        points={{-64,80},{-64,76},{-66,76}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder.frame_a, revolute1.frame_b) annotation (Line(
        points={{-98,112},{-98,106},{-94,106},{-94,100}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder2.frame_a, revolute3.frame_b) annotation (Line(
        points={{-56,112},{-56,106},{-64,106},{-64,100}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute4.frame_a, fixedTranslation3.frame_b) annotation (Line(
        points={{-126,-44},{-114,-44}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder3.frame_a, revolute4.frame_b) annotation (Line(
        points={{-156,-44},{-146,-44}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation6.frame_a, fixedTranslation1.frame_b) annotation (
        Line(
        points={{-102,-18},{-90,-18},{-90,-30},{-74,-30},{-74,-42}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation7.frame_a, fixed.frame_b) annotation (Line(
        points={{-62,-66},{-72,-66},{-72,-68},{-74,-68},{-74,-120},{-128,-120},
            {-128,-138}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute6.frame_a, fixedTranslation1.frame_b) annotation (Line(
        points={{-48,-12},{-58,-12},{-58,-34},{-74,-34},{-74,-42}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute5.frame_a, fixedTranslation1.frame_b) annotation (Line(
        points={{-48,-40},{-60,-40},{-60,-34},{-74,-34},{-74,-42}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation8.frame_a, revolute5.frame_b) annotation (Line(
        points={{-10,-30},{-18,-30},{-18,-38},{-28,-38},{-28,-40}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute5.frame_b, fixedTranslation9.frame_a) annotation (Line(
        points={{-28,-40},{-28,-60},{-14,-60}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation11.frame_a, revolute6.frame_b) annotation (Line(
        points={{-8,-4},{-22,-4},{-22,-12},{-28,-12}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation10.frame_a, revolute6.frame_b) annotation (Line(
        points={{-4,26},{-16,26},{-16,22},{-28,22},{-28,-12}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation11.frame_b, revolute7.frame_a) annotation (Line(
        points={{12,-4},{24,-4}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation9.frame_b, revolute8.frame_a) annotation (Line(
        points={{6,-60},{18,-60}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation12.frame_a, revolute7.frame_b) annotation (Line(
        points={{72,14},{50,14},{50,-4},{44,-4}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute7.frame_b, fixedTranslation13.frame_a) annotation (Line(
        points={{44,-4},{66,-4},{66,-10},{72,-10}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation14.frame_a, revolute8.frame_b) annotation (Line(
        points={{60,-36},{44,-36},{44,-60},{38,-60}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute8.frame_b, fixedTranslation15.frame_a) annotation (Line(
        points={{38,-60},{38,-80},{62,-80},{62,-88}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation12.frame_b, bodyCylinder4.frame_a) annotation (Line(
        points={{92,14},{108,14}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation13.frame_b, bodyCylinder5.frame_a) annotation (Line(
        points={{92,-10},{92,-8},{108,-8}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation14.frame_b, bodyCylinder6.frame_a) annotation (Line(
        points={{80,-36},{80,-34},{98,-34}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder7.frame_a, fixedTranslation15.frame_b) annotation (Line(
        points={{106,-72},{86,-72},{86,-88},{82,-88}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation.frame_a, fixedTranslation16.frame_a) annotation (
        Line(
        points={{-50,-124},{-50,-120},{-14,-120},{-14,-128}},
        color={95,95,95},
        thickness=0.5));
    connect(springDamperParallel.frame_a, fixedTranslation16.frame_b)
      annotation (Line(
        points={{4,-112},{6,-112},{6,-128}},
        color={95,95,95},
        thickness=0.5));
    connect(springDamperParallel.frame_b, fixedTranslation7.frame_b)
      annotation (Line(
        points={{4,-92},{4,-74},{-36,-74},{-36,-66},{-42,-66}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation8.frame_b, springDamperParallel1.frame_a)
      annotation (Line(
        points={{10,-30},{26,-30},{26,36},{48,36},{48,50}},
        color={95,95,95},
        thickness=0.5));
    connect(springDamperParallel1.frame_b, fixedTranslation10.frame_b)
      annotation (Line(
        points={{48,70},{36,70},{36,66},{16,66},{16,26}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation6.frame_b, springDamperParallel2.frame_a)
      annotation (Line(
        points={{-122,-18},{-138,-18},{-138,-4},{-156,-4}},
        color={95,95,95},
        thickness=0.5));
    connect(springDamperParallel2.frame_b, fixedTranslation17.frame_b)
      annotation (Line(
        points={{-156,16},{-152,16},{-152,26},{-120,26},{-120,16}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute9.frame_a, fixedTranslation2.frame_b) annotation (Line(
        points={{-87,-4},{-82,-4},{-82,-8},{-74,-8}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation17.frame_a, revolute9.frame_b) annotation (Line(
        points={{-120,36},{-112,36},{-112,10},{-87,10}},
        color={95,95,95},
        thickness=0.5));
    connect(prismatic.frame_a, revolute9.frame_b) annotation (Line(
        points={{-86,20},{-86,15},{-87,15},{-87,10}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation4.frame_a, prismatic.frame_b) annotation (Line(
        points={{-96,56},{-96,40},{-86,40}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation5.frame_a, prismatic.frame_b) annotation (Line(
        points={{-66,56},{-66,46},{-86,46},{-86,40}},
        color={95,95,95},
        thickness=0.5));
    connect(position.flange, prismatic.axis) annotation (Line(points={{-128,66},
            {-124,66},{-124,50},{-92,50},{-92,38}}, color={0,127,0}));
    connect(realExpression.y, position.s_ref)
      annotation (Line(points={{-171,66},{-150,66}}, color={0,0,127}));
    connect(contact1.frame_b, fixedTranslation12.frame_b) annotation (Line(
        points={{154.2,30.4},{92,30.4},{92,14}},
        color={95,95,95},
        thickness=0.5));
    connect(contact2.frame_b, fixedTranslation13.frame_b) annotation (Line(
        points={{154.2,-7.6},{144,-7.6},{144,-22},{92,-22},{92,-10}},
        color={95,95,95},
        thickness=0.5));
    connect(contact1.frame_a, cinghia1.frame_b) annotation (Line(
        points={{176,30.2},{194,30.2},{194,18},{268,18},{268,-22.52}},
        color={95,95,95},
        thickness=0.5));
    connect(contact2.frame_a, cinghia1.frame_b) annotation (Line(
        points={{176,-7.8},{192,-7.8},{192,-12},{268,-12},{268,-22.52}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation14.frame_b, contact3.frame_b) annotation (Line(
        points={{80,-36},{82,-36},{82,-58},{150.2,-58},{150.2,-39.6}},
        color={95,95,95},
        thickness=0.5));
    connect(contact3.frame_a, cinghia1.frame_b) annotation (Line(
        points={{172,-39.8},{204,-39.8},{204,-34},{268,-34},{268,-22.52}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation15.frame_b, contact4.frame_b) annotation (Line(
        points={{82,-88},{100,-88},{100,-100},{142.2,-100},{142.2,-81.6}},
        color={95,95,95},
        thickness=0.5));
    connect(contact4.frame_a, cinghia1.frame_b) annotation (Line(
        points={{164,-81.8},{208,-81.8},{208,-68},{268,-68},{268,-22.52}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute3.frame_b, contact5.frame_b) annotation (Line(
        points={{-64,100},{-46,100},{-46,104},{-13.8,104},{-13.8,106.4}},
        color={95,95,95},
        thickness=0.5));
    connect(contact5.frame_a, cinghia1.frame_b) annotation (Line(
        points={{8,106.2},{118,106.2},{118,74},{268,74},{268,-22.52}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute1.frame_b, contact6.frame_b) annotation (Line(
        points={{-94,100},{-94,162.4},{-45.8,162.4}},
        color={95,95,95},
        thickness=0.5));
    connect(contact6.frame_a, cinghia1.frame_b) annotation (Line(
        points={{-24,162.2},{44,162.2},{44,142},{268,142},{268,-22.52}},
        color={95,95,95},
        thickness=0.5));
    connect(contact7.frame_b, revolute4.frame_b) annotation (Line(
        points={{-190.2,-29.6},{-146,-29.6},{-146,-44}},
        color={95,95,95},
        thickness=0.5));
    connect(contact7.frame_a, cinghia1.frame_b) annotation (Line(
        points={{-212,-29.8},{-252,-29.8},{-252,32},{-236,32},{-236,194},{268,
            194},{268,-22.52}},
        color={95,95,95},
        thickness=0.5));
    connect(contact_sprocket1.frame_a, cinghia1.frame_b) annotation (Line(
        points={{122,-185.8},{178,-185.8},{178,-144},{268,-144},{268,-22.52}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute2.frame_b, contact_sprocket1.frame_b) annotation (Line(
        points={{-62,-179},{-2,-179},{-2,-216},{100.2,-216},{100.2,-185.6}},
        color={95,95,95},
        thickness=0.5));
    connect(rotational_stop1.frame_a, revolute7.frame_b) annotation (Line(
        points={{46,19},{50,19},{50,-4},{44,-4}},
        color={95,95,95},
        thickness=0.5));
    connect(rotational_stop1.frame_b, revolute7.frame_a) annotation (Line(
        points={{30,19},{24,19},{24,-4}},
        color={95,95,95},
        thickness=0.5));
    connect(rotational_stop2.frame_a, revolute8.frame_b) annotation (Line(
        points={{40,-37},{42,-37},{42,-60},{38,-60}},
        color={95,95,95},
        thickness=0.5));
    connect(rotational_stop2.frame_b, revolute8.frame_a) annotation (Line(
        points={{24,-37},{20,-37},{20,-38},{18,-38},{18,-60}},
        color={95,95,95},
        thickness=0.5));
    connect(position1.flange, revolute.axis) annotation (Line(points={{-136,-76},
            {-68,-76},{-68,-70},{-54,-70}}, color={0,0,0}));
    connect(realExpression1.y, position1.phi_ref) annotation (Line(points={{
            -189,-76},{-180,-76},{-180,-78},{-158,-78},{-158,-76}}, color={0,0,
            127}));
    connect(realExpression2.y, position2.phi_ref) annotation (Line(points={{
            -221,44},{-212,44},{-212,42},{-190,42},{-190,44}}, color={0,0,127}));
    connect(position2.flange, revolute9.axis) annotation (Line(points={{-168,44},
            {-156,44},{-156,38},{-126,38},{-126,3},{-94,3}}, color={0,0,0}));
    connect(realExpression3.y, position3.phi_ref) annotation (Line(points={{-37,
            42},{-32,42},{-32,20},{-44,20},{-44,10},{-70,10},{-70,24}}, color={
            0,0,127}));
    connect(position3.flange, revolute6.axis) annotation (Line(points={{-48,24},
            {-42,24},{-42,14},{-38,14},{-38,-2}}, color={0,0,0}));
    connect(realExpression4.y, position4.phi_ref) annotation (Line(points={{
            -199,-106},{-120,-106},{-120,-100},{-148,-100}}, color={0,0,127}));
    connect(position4.flange, revolute5.axis) annotation (Line(points={{-126,
            -100},{-114,-100},{-114,-98},{-38,-98},{-38,-30}}, color={0,0,0}));
    connect(realExpression5.y, position5.phi_ref) annotation (Line(points={{63,
            -126},{72,-126},{72,-128},{94,-128},{94,-126}}, color={0,0,127}));
    connect(position5.flange, revolute8.axis) annotation (Line(points={{116,
            -126},{132,-126},{132,-120},{28,-120},{28,-50}}, color={0,0,0}));
    connect(realExpression6.y, position6.phi_ref) annotation (Line(points={{87,
            54},{96,54},{96,52},{118,52},{118,54}}, color={0,0,127}));
    connect(position6.flange, revolute7.axis) annotation (Line(points={{140,54},
            {148,54},{148,6},{34,6}}, color={0,0,0}));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=30,
        Tolerance=0.001,
        __Dymola_Algorithm="Dassl"));
  end track_bunker_pro;

  model rotational_stop
    Modelica.Units.SI.Angle phi;
    Modelica.Units.SI.AngularVelocity w;
    parameter Real K = 1e8;
    parameter Real D = 1e5;
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b
      annotation (Placement(transformation(extent={{-116,-16},{-84,16}})));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a
      annotation (Placement(transformation(extent={{84,-16},{116,16}})));
  equation
    w = frame_a.R.w[3] - frame_b.R.w[3];
    der(phi) = w;

    if noEvent(phi > (3/180)*Modelica.Constants.pi) then
      frame_a.t[3] = max(0, K*(phi - (3/180)*Modelica.Constants.pi) + D*w);
    elseif noEvent(phi < -(3/180)*Modelica.Constants.pi) then
      frame_a.t[3] = min(0, K*(phi + (3/180)*Modelica.Constants.pi) + D*w);
    else
      frame_a.t[3] = 0;
    end if;
    frame_a.t[1] = 0;
    frame_a.t[2] = 0;
    frame_a.f = {0,0,0};
    frame_a.f + frame_b.f = {0, 0, 0};

   frame_a.t + frame_b.t = {0, 0, 0};

    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end rotational_stop;

  model trackterrain_grouser
    parameter Boolean anim = true;
    parameter Modelica.Units.SI.Position start[3] = {0, 0, 0};
     parameter Modelica.Units.SI.Angle a = 0;
    parameter Modelica.Units.SI.Angle b = 0;
    parameter Modelica.Units.SI.Angle g = 0;
    final parameter Modelica.Mechanics.MultiBody.Frames.Orientation R = Modelica.Mechanics.MultiBody.Frames.axesRotations({1,2,3}, {a,b,g}, zeros(3));
    final parameter Modelica.Units.SI.Position start1[3] = start + Modelica.Mechanics.MultiBody.Frames.resolve1(R, - {-0.18554,-0.07342,0} + {0.21488,-0.04227,0} + {0.22252,0.01926,0} + {0.20314,0.02838,0});
    input Modelica.Units.SI.Length dh[80];
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute(useAxisFlange=false,
        animation=false,
      phi(fixed=true, displayUnit="rad"))
      annotation (Placement(transformation(extent={{-64,-90},{-44,-70}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation(animation=
         anim,                                                           r={-0.18554,
          -0.07342,0},
      color={0,0,0})     annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-52,-128})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute2(useAxisFlange=true,
        animation=false) annotation (Placement(transformation(
          extent={{-10,11},{10,-11}},
          rotation=270,
          origin={-92,-143})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder Sprocket(
      animation=anim,
      r={0,0,0},
      r_shape={0,0,-0.105*0.5},
      length=0.105,
      diameter=0.16,
      density=1399.4110125,
      color={0,0,0},
      r_0(start=start),
      angles_fixed=false,
      angles_start={a,b,g},
      useQuaternions=false) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-44,-238})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation1(animation=
         anim,                                                            r={
          0.21488,-0.04227,0},
      color={0,0,0})       annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-74,-52})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation2(animation=
         anim,                                                            r={
          0.22252,0.01926,0},
      color={0,0,0})        annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-74,-18})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation3(animation=
         anim,                                                            r={-0.03379,
          0.15530,0},
      color={0,0,0})
      annotation (Placement(transformation(extent={{-94,-54},{-114,-34}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation4(animation=
         anim,                                                            r={
          0.20314,0.02838,0},
      color={0,0,0})      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-96,66})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation5(animation=
         anim,                                                            r={
          0.10492,-0.09153,0},
      color={0,0,0})       annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-66,66})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute1(animation=false)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-150,84})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute3(animation=false)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-70,92})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder Idler_1(
      animation=anim,
      r={0,0,0},
      r_shape={0,0,-0.15*0.5},
      length=0.15,
      diameter=0.139,
      density(displayUnit="kg/m3") = 722.8451943481186,
      color={0,0,0},
      r_0(start=start1),
      angles_fixed=false,
      angles_start={a,b,g},
      useQuaternions=false) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-150,134})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder Idler_2(
      animation=anim,
      r={0,0,0},
      r_shape={0,0,-0.15*0.5},
      length=0.15,
      diameter=0.126,
      density(displayUnit="kg/m3") = 794.5646258503402,
      color={0,0,0},
      angles_fixed=false,
      angles_start={a,b,g},
      useQuaternions=false) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-82,134})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute4(animation=false)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={-134,-78})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder Upper_wheel(
      animation=anim,
      r={0,0,0},
      r_shape={0,0,-0.15*0.5},
      length=0.15,
      diameter=0.08,
      density=1142.9184375,
      color={0,0,0},
      angles_fixed=false,
      angles_start={a,b,g},
      useQuaternions=false)
      annotation (Placement(transformation(extent={{-182,-80},{-202,-60}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation6(animation=
         anim,                                                            r={
          0.17068,0.08769,0},
      color={0,0,0})
      annotation (Placement(transformation(extent={{-102,-28},{-122,-8}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation7(animation=
         anim,                                                            r={
          0.06426,-0.07476,0},
      color={0,0,0})
      annotation (Placement(transformation(extent={{-62,-76},{-42,-56}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute5(useAxisFlange=false,
        animation=false,
      phi(fixed=true))
      annotation (Placement(transformation(extent={{-48,-50},{-28,-30}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute6(useAxisFlange=false,
        animation=false,
      phi(fixed=true))
      annotation (Placement(transformation(extent={{-48,-22},{-28,-2}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation8(animation=
         anim,                                                            r={
          0.07501,0.07623,0},
      color={0,0,0})
      annotation (Placement(transformation(extent={{-10,-40},{10,-20}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation9(animation=
         anim,                                                            r={
          0.11355,-0.09018,0},
      color={0,0,0})
      annotation (Placement(transformation(extent={{-14,-70},{6,-50}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation10(animation=
         anim,                                                             r={-0.07501,
          0.07623,0},
      color={0,0,0})
      annotation (Placement(transformation(extent={{-24,44},{-4,64}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation11(animation=
         anim,                                                             r={-0.11355,
          -0.09018,0},
      color={0,0,0})
      annotation (Placement(transformation(extent={{-8,-14},{12,6}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute7(useAxisFlange=false,
        animation=false,
      phi(fixed=true))
      annotation (Placement(transformation(extent={{38,-14},{58,6}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute8(useAxisFlange=false,
        animation=false,
      phi(fixed=true))
      annotation (Placement(transformation(extent={{34,-90},{54,-70}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation12(animation=
         anim,                                                             r={
          0.05250,-0.01633,0},
      color={0,0,0})
      annotation (Placement(transformation(extent={{74,48},{94,68}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation13(animation=
         anim,                                                             r={-0.05250,
          -0.01633,0},
      color={0,0,0})
      annotation (Placement(transformation(extent={{80,-16},{100,4}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation14(animation=
         anim,                                                             r={
          0.05250,-0.01633,0},
      color={0,0,0})
      annotation (Placement(transformation(extent={{72,-56},{92,-36}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation15(animation=
         anim,                                                             r={-0.05250,
          -0.01633,0},
      color={0,0,0})
      annotation (Placement(transformation(extent={{70,-122},{90,-102}})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder Road_wheel_4(
      animation=anim,
      r={0,0,0},
      r_shape={0,0,-0.024*0.5},
      length=0.024,
      diameter=0.069,
      density=2295.0850661626,
      color={0,0,0},
      angles_fixed=false,
      angles_start={a,b,g},
      w_0_fixed=false,
      useQuaternions=false)
      annotation (Placement(transformation(extent={{152,48},{172,68}})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder Road_wheel_3(
      animation=anim,
      r={0,0,0},
      r_shape={0,0,-0.024*0.5},
      length=0.024,
      diameter=0.069,
      density=2295.0850661626,
      color={0,0,0},
      angles_fixed=false,
      angles_start={a,b,g},
      w_0_fixed=false,
      useQuaternions=false)
      annotation (Placement(transformation(extent={{154,-16},{174,4}})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder Road_wheel_2(
      animation=anim,
      r={0,0,0},
      r_shape={0,0,-0.024*0.5},
      length=0.024,
      diameter=0.069,
      density=2295.0850661626,
      color={0,0,0},
      angles_fixed=false,
      angles_start={a,b,g},
      w_0_fixed=false,
      useQuaternions=false)
      annotation (Placement(transformation(extent={{154,-54},{174,-34}})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder Road_wheel_1(
      animation=anim,
      r={0,0,0},
      r_shape={0,0,-0.024*0.5},
      length=0.024,
      diameter=0.069,
      density=2295.0850661626,
      color={0,0,0},
      angles_fixed=false,
      angles_start={a,b,g},
      w_0_fixed=false,
      useQuaternions=false)
      annotation (Placement(transformation(extent={{154,-118},{174,-98}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation16(animation=
         anim,                                                             r={
          0.084035,0.0331024,0},
      color={0,0,0})
      annotation (Placement(transformation(extent={{-14,-138},{6,-118}})));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel
      springDamperParallel(
      animation=anim,
      c=131389.9743,
      s_unstretched=0.10,
      d=13138.99743)
              annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={4,-102})));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel
      springDamperParallel1(
      animation=anim,
      c=131389.9743,
      s_unstretched=0.15375,
      d=13138.99743)
              annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={28,72})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation17(animation=
         anim,                                                             r={
          0.20314,0.02838,0} + {-0.12,0.04005,0},
      color={0,0,0})                       annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-130,6})));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel
      springDamperParallel2(
      animation=anim,
      c=131389.9743,
      s_unstretched=0.134,
      d=13138.99743)
              annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-170,-14})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute9(useAxisFlange=false,
        animation=false,
      phi(fixed=true))   annotation (Placement(transformation(
          extent={{-7,-7},{7,7}},
          rotation=90,
          origin={-87,3})));
    Modelica.Mechanics.MultiBody.Joints.Prismatic prismatic(useAxisFlange=true,
        animation=false)
                       annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-86,30})));
    Modelica.Mechanics.Translational.Sources.Position position(exact=true)
      annotation (Placement(transformation(extent={{-156,32},{-136,52}})));
    Modelica.Blocks.Sources.RealExpression realExpression
      annotation (Placement(transformation(extent={{-230,20},{-210,40}})));
    contact contact1(
      R=0.069*0.5,
      lar=0.024,
      L=0.024)
      annotation (Placement(transformation(extent={{216,44},{236,64}})));
    contact contact2(
      R=0.069*0.5,
      lar=0.024,
      L=0.024)
      annotation (Placement(transformation(extent={{214,-16},{234,4}})));
    contact contact3(
      R=0.069*0.5,
      lar=0.024,
      L=0.024)
      annotation (Placement(transformation(extent={{212,-54},{232,-34}})));
    contact contact4(
      R=0.069*0.5,
      lar=0.024,
      L=0.024)
      annotation (Placement(transformation(extent={{216,-114},{236,-94}})));
    contact contact5(R=0.126*0.5)
      annotation (Placement(transformation(extent={{-4,142},{16,162}})));
    contact contact6(R=0.139*0.5)
      annotation (Placement(transformation(extent={{-42,174},{-22,194}})));
    contact contact7(R=0.08*0.5)
      annotation (Placement(transformation(extent={{-236,-64},{-256,-44}})));
    contact_sprocket contact_sprocket1
      annotation (Placement(transformation(extent={{102,-196},{122,-176}})));
    rotational_stop rotational_stop1(K=1e6, D=10)
      annotation (Placement(transformation(extent={{42,16},{58,30}})));
    rotational_stop rotational_stop2(K=1e6, D=10)
      annotation (Placement(transformation(extent={{34,-56},{50,-42}})));
     Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a annotation (
      Placement(visible = true, transformation(origin={-310,-158},   extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-100, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b[80] annotation (
      Placement(visible = true, transformation(origin={482,-8},    extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {98, 0}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_a annotation (
        Placement(transformation(extent={{-10,202},{10,222}}),
                                                             iconTransformation(
            extent={{-10,202},{10,222}})));
    track_with_grouser track(
      anim=anim,
      start=start,
      dh=dh,
      l=0.0586*38,
      startX={0.13524441,0.16336884,0.19149324,0.21961766,0.24774207,0.27586648,
          0.3039909,0.33211532,0.3602397,0.38836285,0.41647583,0.44458714,
          0.47269815,0.5008085,0.52891237,0.55688405,0.58262175,0.6032808,
          0.62217665,0.6408851,0.65957004,0.6782187,0.696536,0.7112605,0.716824,
          0.71197385,0.69739616,0.67507076,0.64832556,0.62066644,0.59294254,
          0.56521285,0.537483,0.5097535,0.48202458,0.4542961,0.42656812,
          0.39884064,0.37111363,0.3433871,0.31566107,0.28793547,0.26020995,
          0.23248067,0.20471112,0.17665297,0.1495237,0.123243846,0.09706277,
          0.07089086,0.04471898,0.018546196,-0.007627586,-0.03380237,-0.059978165,
          -0.08615496,-0.11233275,-0.13851139,-0.1646895,-0.19085327,-0.21686272,
          -0.24097931,-0.25945207,-0.27013397,-0.2719339,-0.26467142,-0.24910323,
          -0.226839,-0.20030162,-0.17248595,-0.14458823,-0.11668263,-0.0887759,
          -0.06086873,-0.032961085,-0.0050521013,0.022866843,0.050872553,
          0.07899562,0.10712001} + fill(0.18554, 80),
      startY={-0.19040328,-0.19041218,-0.19041839,-0.19042167,-0.19042195,-0.19041924,
          -0.19041334,-0.19040224,-0.19036625,-0.19010034,-0.18929689,-0.18843658,
          -0.1875661,-0.18667395,-0.18559812,-0.18267798,-0.17135578,-0.15227471,
          -0.1314442,-0.11044516,-0.08942516,-0.068372965,-0.04703182,-0.023075268,
          0.004486749,0.032183215,0.05622712,0.073321044,0.08200841,0.08709987,
          0.09182613,0.096518226,0.101209484,0.10590327,0.11059991,0.11529946,
          0.120001905,0.12470725,0.1294155,0.13412665,0.13884066,0.14355734,
          0.1482744,0.1529695,0.15742008,0.15932229,0.15191655,0.14190124,
          0.13163047,0.12133636,0.11104223,0.10075044,0.0904612,0.08017456,
          0.06989051,0.05960904,0.049330126,0.03905341,0.02877536,0.018460942,
          0.007763456,-0.006696874,-0.027895795,-0.05390605,-0.08196663,-0.10913077,
          -0.1325458,-0.14971939,-0.15901689,-0.16316524,-0.16672072,-0.17021371,
          -0.17369777,-0.17717819,-0.18065481,-0.1841206,-0.18750517,-0.19007342,
          -0.19034702,-0.19038792} + fill(0.07342, 80),
      startZ=fill(0, 80),
      thZ(displayUnit="rad") = {-0.0003167978,-0.00022154418,-0.00011671506,-1.08526265e-05,
        9.582425e-05,0.00020965889,0.00039465772,0.0012881175,0.009234589,
        0.028549504,0.03059104,0.030956108,0.0317336,0.038335502,0.10310366,
        0.4122258,0.7446929,0.8339565,0.8430107,0.8441416,0.8458825,0.861678,
        1.0180945,1.3697258,1.7422643,2.1139288,2.4861243,2.8260055,2.9593987,
        2.972725,2.9739714,2.9740026,2.9739118,2.9738085,2.973704,2.9735994,
        2.973495,2.9733903,2.973286,2.9731822,2.9730868,2.9730735,2.973875,
        2.9827793,3.0722399,-2.8762295,-2.7775958,-2.7677548,-2.766853,-2.7668512,
        -2.7669399,-2.7670376,-2.7671363,-2.767235,-2.7673337,-2.7674313,-2.7675147,
        -2.7674623,-2.766057,-2.7512207,-2.6029167,-2.2892492,-1.9621332,-1.6365067,
        -1.3112022,-0.9856874,-0.65874875,-0.33877757,-0.14829077,-0.12678835,-0.12452726,
        -0.124204345,-0.12407412,-0.1239377,-0.12354662,-0.120608404,-0.092386425,
        -0.00982342,-0.0014643206,-0.00054822356},
      there_is_grouser={true,true,true,true,true,false,true,false,true,false,
          true,false,true,false,true,false,true,true,true,true,true,true,true,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,true,true,true,true,true,true,true,true,true,false,true,false,
          true,false,true,false,true,false,true,false,true,true,true,true,true,
          true,true,true,false,true,false,true,false,true,false,true,false,true,
          false,true,true,true,true,true},
      cm={0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,0.01172,0,
          0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,0.0029175,
          0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,0.010255,0,
          0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,0.002185,
          0.024185,0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,
          0.01172,0,0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,
          0.0029175,0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,
          0.010255,0,0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,
          0.002185,0.024185},
      lg={0.0073,0.00437,0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,
          0.0146,0,0.0146,0,0.011695,0.002905,0.008765,0.005835,0.005835,
          0.008765,0.002905,0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,
          0.0146,0,0.01316,0.00144,0.01023,0.00437,0.0073,0.0073,0.00437,
          0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,
          0,0.011695,0.002905,0.008765,0.005835,0.005835,0.008765,0.002905,
          0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.01316,
          0.00144,0.01023,0.00437,0.0073},
      a=a,
      b=b,
      g=g) annotation (Placement(visible=true, transformation(
          origin={367,-3},
          extent={{-67,-67},{67,67}},
          rotation=0)));

    chain chain1
      annotation (Placement(transformation(extent={{-84,-200},{-64,-180}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute10(useAxisFlange=false,
        animation=false)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-34,-184})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute11(useAxisFlange=false,
        animation=false)
      annotation (Placement(transformation(extent={{108,48},{128,68}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute12(useAxisFlange=false,
        animation=false)
      annotation (Placement(transformation(extent={{110,-16},{130,4}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute13(useAxisFlange=false,
        animation=false)
      annotation (Placement(transformation(extent={{112,-52},{132,-32}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute14(useAxisFlange=false,
        animation=false)
      annotation (Placement(transformation(extent={{112,-118},{132,-98}})));
  equation
    connect(revolute.frame_b, fixedTranslation.frame_a) annotation (Line(
        points={{-44,-80},{-38,-80},{-38,-118},{-52,-118}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation2.frame_a, fixedTranslation1.frame_b) annotation (
        Line(
        points={{-74,-28},{-74,-42}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation3.frame_a, fixedTranslation1.frame_b) annotation (
        Line(
        points={{-94,-44},{-88,-44},{-88,-34},{-74,-34},{-74,-42}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute1.frame_a, fixedTranslation4.frame_b) annotation (Line(
        points={{-150,74},{-112,74},{-112,76},{-96,76}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute3.frame_a, fixedTranslation5.frame_b) annotation (Line(
        points={{-70,82},{-70,76},{-66,76}},
        color={95,95,95},
        thickness=0.5));
    connect(Idler_1.frame_a, revolute1.frame_b) annotation (Line(
        points={{-150,124},{-150,94}},
        color={95,95,95},
        thickness=0.5));
    connect(Idler_2.frame_a, revolute3.frame_b) annotation (Line(
        points={{-82,124},{-82,118},{-70,118},{-70,102}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute4.frame_a, fixedTranslation3.frame_b) annotation (Line(
        points={{-124,-78},{-120,-78},{-120,-44},{-114,-44}},
        color={95,95,95},
        thickness=0.5));
    connect(Upper_wheel.frame_a, revolute4.frame_b) annotation (Line(
        points={{-182,-70},{-160,-70},{-160,-78},{-144,-78}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation6.frame_a, fixedTranslation1.frame_b) annotation (
        Line(
        points={{-102,-18},{-90,-18},{-90,-30},{-74,-30},{-74,-42}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute6.frame_a, fixedTranslation1.frame_b) annotation (Line(
        points={{-48,-12},{-58,-12},{-58,-34},{-74,-34},{-74,-42}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute5.frame_a, fixedTranslation1.frame_b) annotation (Line(
        points={{-48,-40},{-60,-40},{-60,-34},{-74,-34},{-74,-42}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation8.frame_a, revolute5.frame_b) annotation (Line(
        points={{-10,-30},{-18,-30},{-18,-38},{-28,-38},{-28,-40}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute5.frame_b, fixedTranslation9.frame_a) annotation (Line(
        points={{-28,-40},{-28,-60},{-14,-60}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation11.frame_a, revolute6.frame_b) annotation (Line(
        points={{-8,-4},{-22,-4},{-22,-12},{-28,-12}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation10.frame_a, revolute6.frame_b) annotation (Line(
        points={{-24,54},{-28,54},{-28,2},{-24,2},{-24,-6},{-22,-6},{-22,-12},{
            -28,-12}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation11.frame_b, revolute7.frame_a) annotation (Line(
        points={{12,-4},{38,-4}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation9.frame_b, revolute8.frame_a) annotation (Line(
        points={{6,-60},{28,-60},{28,-80},{34,-80}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation12.frame_a, revolute7.frame_b) annotation (Line(
        points={{74,58},{72,58},{72,-4},{58,-4}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute7.frame_b, fixedTranslation13.frame_a) annotation (Line(
        points={{58,-4},{72,-4},{72,-8},{76,-8},{76,-6},{80,-6}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation14.frame_a, revolute8.frame_b) annotation (Line(
        points={{72,-46},{58,-46},{58,-80},{54,-80}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute8.frame_b, fixedTranslation15.frame_a) annotation (Line(
        points={{54,-80},{64,-80},{64,-112},{70,-112}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation.frame_a, fixedTranslation16.frame_a) annotation (
        Line(
        points={{-52,-118},{-52,-120},{-14,-120},{-14,-128}},
        color={95,95,95},
        thickness=0.5));
    connect(springDamperParallel.frame_a, fixedTranslation16.frame_b)
      annotation (Line(
        points={{4,-112},{6,-112},{6,-128}},
        color={95,95,95},
        thickness=0.5));
    connect(springDamperParallel.frame_b, fixedTranslation7.frame_b)
      annotation (Line(
        points={{4,-92},{4,-74},{-36,-74},{-36,-66},{-42,-66}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation8.frame_b, springDamperParallel1.frame_a)
      annotation (Line(
        points={{10,-30},{16,-30},{16,54},{28,54},{28,62}},
        color={95,95,95},
        thickness=0.5));
    connect(springDamperParallel1.frame_b, fixedTranslation10.frame_b)
      annotation (Line(
        points={{28,82},{28,84},{2,84},{2,54},{-4,54}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation6.frame_b, springDamperParallel2.frame_a)
      annotation (Line(
        points={{-122,-18},{-138,-18},{-138,-24},{-170,-24}},
        color={95,95,95},
        thickness=0.5));
    connect(springDamperParallel2.frame_b, fixedTranslation17.frame_b)
      annotation (Line(
        points={{-170,-4},{-150,-4},{-150,-4},{-130,-4}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute9.frame_a, fixedTranslation2.frame_b) annotation (Line(
        points={{-87,-4},{-82,-4},{-82,-8},{-74,-8}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation17.frame_a, revolute9.frame_b) annotation (Line(
        points={{-130,16},{-90,16},{-90,10},{-87,10}},
        color={95,95,95},
        thickness=0.5));
    connect(prismatic.frame_a, revolute9.frame_b) annotation (Line(
        points={{-86,20},{-86,15},{-87,15},{-87,10}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation4.frame_a, prismatic.frame_b) annotation (Line(
        points={{-96,56},{-96,40},{-86,40}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation5.frame_a, prismatic.frame_b) annotation (Line(
        points={{-66,56},{-66,46},{-86,46},{-86,40}},
        color={95,95,95},
        thickness=0.5));
    connect(position.flange, prismatic.axis) annotation (Line(points={{-136,42},
            {-120,42},{-120,38},{-92,38}},          color={0,127,0}));
    connect(realExpression.y, position.s_ref)
      annotation (Line(points={{-209,30},{-160,30},{-160,42},{-158,42}},
                                                     color={0,0,127}));
    connect(rotational_stop1.frame_a, revolute7.frame_b) annotation (Line(
        points={{58,23},{60,23},{60,-4},{58,-4}},
        color={95,95,95},
        thickness=0.5));
    connect(rotational_stop1.frame_b, revolute7.frame_a) annotation (Line(
        points={{42,23},{42,24},{28,24},{28,-4},{38,-4}},
        color={95,95,95},
        thickness=0.5));
    connect(rotational_stop2.frame_a, revolute8.frame_b) annotation (Line(
        points={{50,-49},{58,-49},{58,-80},{54,-80}},
        color={95,95,95},
        thickness=0.5));
    connect(rotational_stop2.frame_b, revolute8.frame_a) annotation (Line(
        points={{34,-49},{30,-49},{30,-74},{28,-74},{28,-80},{34,-80}},
        color={95,95,95},
        thickness=0.5));

    connect(contact1.frame_a, track.frame_b) annotation (Line(
        points={{236,54.2},{274,54.2},{274,56},{300,56},{300,-4.34}},
        color={95,95,95},
        thickness=0.5));
    connect(contact2.frame_a, track.frame_b) annotation (Line(
        points={{234,-5.8},{267,-5.8},{267,-4.34},{300,-4.34}},
        color={95,95,95},
        thickness=0.5));
    connect(contact3.frame_a, track.frame_b) annotation (Line(
        points={{232,-43.8},{270,-43.8},{270,-28},{300,-28},{300,-4.34}},
        color={95,95,95},
        thickness=0.5));
    connect(contact4.frame_a, track.frame_b) annotation (Line(
        points={{236,-103.8},{248,-103.8},{248,-64},{300,-64},{300,-4.34}},
        color={95,95,95},
        thickness=0.5));
    connect(contact_sprocket1.frame_a, track.frame_b) annotation (Line(
        points={{122,-185.8},{162,-185.8},{162,-164},{300,-164},{300,-4.34}},
        color={95,95,95},
        thickness=0.5));
    connect(contact5.frame_a, track.frame_b) annotation (Line(
        points={{16,152.2},{74,152.2},{74,102},{300,102},{300,-4.34}},
        color={95,95,95},
        thickness=0.5));
    connect(contact6.frame_a, track.frame_b) annotation (Line(
        points={{-22,184.2},{106,184.2},{106,138},{300,138},{300,-4.34}},
        color={95,95,95},
        thickness=0.5));
    connect(contact7.frame_a, track.frame_b) annotation (Line(
        points={{-256,-53.8},{-256,204},{300,204},{300,-4.34}},
        color={95,95,95},
        thickness=0.5));
    connect(track.frame_a, frame_b) annotation (Line(
        points={{434,-8.36},{458,-8.36},{458,-8},{482,-8}},
        color={95,95,95},
        thickness=0.5));
    connect(flange_a, revolute2.axis) annotation (Line(points={{0,212},{-274,
            212},{-274,-202},{-103,-202},{-103,-143}},     color={0,0,0}));
    connect(revolute.frame_a, revolute2.frame_a) annotation (Line(
        points={{-64,-80},{-76,-80},{-76,-118},{-92,-118},{-92,-133}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute2.frame_b, chain1.frame_a) annotation (Line(
        points={{-92,-153},{-92,-190},{-84.2,-190}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation.frame_b, revolute10.frame_a) annotation (Line(
        points={{-52,-138},{-46,-138},{-46,-164},{-34,-164},{-34,-174}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute10.frame_b, Sprocket.frame_a) annotation (Line(
        points={{-34,-194},{-34,-214},{-44,-214},{-44,-228}},
        color={95,95,95},
        thickness=0.5));
    connect(chain1.frame_b, Sprocket.frame_b) annotation (Line(
        points={{-64.2,-190},{-60,-190},{-60,-250},{-44,-250},{-44,-248}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation12.frame_b, revolute11.frame_a) annotation (Line(
        points={{94,58},{108,58}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute11.frame_b, Road_wheel_4.frame_a) annotation (Line(
        points={{128,58},{152,58}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation13.frame_b, revolute12.frame_a) annotation (Line(
        points={{100,-6},{110,-6}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute12.frame_b, Road_wheel_3.frame_a) annotation (Line(
        points={{130,-6},{154,-6}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation14.frame_b, revolute13.frame_a) annotation (Line(
        points={{92,-46},{106,-46},{106,-42},{112,-42}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute13.frame_b, Road_wheel_2.frame_a) annotation (Line(
        points={{132,-42},{143,-42},{143,-44},{154,-44}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation15.frame_b, revolute14.frame_a) annotation (Line(
        points={{90,-112},{90,-108},{112,-108}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute14.frame_b, Road_wheel_1.frame_a) annotation (Line(
        points={{132,-108},{154,-108}},
        color={95,95,95},
        thickness=0.5));
    connect(contact1.frame_b, Road_wheel_4.frame_b) annotation (Line(
        points={{214.2,54.4},{193,54.4},{193,58},{172,58}},
        color={95,95,95},
        thickness=0.5));
    connect(contact2.frame_b, Road_wheel_3.frame_b) annotation (Line(
        points={{212.2,-5.6},{193.1,-5.6},{193.1,-6},{174,-6}},
        color={95,95,95},
        thickness=0.5));
    connect(contact3.frame_b, Road_wheel_2.frame_b) annotation (Line(
        points={{210.2,-43.6},{192.1,-43.6},{192.1,-44},{174,-44}},
        color={95,95,95},
        thickness=0.5));
    connect(Road_wheel_1.frame_b, contact4.frame_b) annotation (Line(
        points={{174,-108},{192,-108},{192,-106},{202,-106},{202,-103.6},{214.2,
            -103.6}},
        color={95,95,95},
        thickness=0.5));
    connect(contact6.frame_b, Idler_1.frame_b) annotation (Line(
        points={{-43.8,184.4},{-43.8,186},{-150,186},{-150,144}},
        color={95,95,95},
        thickness=0.5));
    connect(contact5.frame_b, Idler_2.frame_b) annotation (Line(
        points={{-5.8,152.4},{-5.8,154},{-66,154},{-66,156},{-82,156},{-82,144}},
        color={95,95,95},
        thickness=0.5));
    connect(contact7.frame_b, Upper_wheel.frame_b) annotation (Line(
        points={{-234.2,-53.6},{-226,-53.6},{-226,-54},{-218,-54},{-218,-70},{
            -202,-70}},
        color={95,95,95},
        thickness=0.5));
    connect(Sprocket.frame_b, contact_sprocket1.frame_b) annotation (Line(
        points={{-44,-248},{4,-248},{4,-246},{66,-246},{66,-185.6},{100.2,-185.6}},
        color={95,95,95},
        thickness=0.5));

    connect(fixedTranslation1.frame_a, revolute.frame_a) annotation (Line(
        points={{-74,-62},{-76,-62},{-76,-80},{-64,-80}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation7.frame_a, fixedTranslation1.frame_a) annotation (
        Line(
        points={{-62,-66},{-70,-66},{-70,-62},{-74,-62}},
        color={95,95,95},
        thickness=0.5));
    connect(frame_a, fixedTranslation1.frame_b) annotation (Line(
        points={{-310,-158},{-120,-158},{-120,-58},{-90,-58},{-90,-48},{-88,-48},
            {-88,-34},{-74,-34},{-74,-42}},
        color={95,95,95},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(__Dymola_Algorithm="Cvode"));
  end trackterrain_grouser;

  model trackterrainfinal1_grouser_pro
     parameter Real start[3] = {0, 0, 0};
   input Modelica.Units.SI.Length dh[80];
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute(useAxisFlange=false,
        animation=false)
      annotation (Placement(transformation(extent={{-64,-90},{-44,-70}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation(r={-0.1749185,
          -0.0689022,0}) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-50,-134})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute2(useAxisFlange=true,
        animation=false) annotation (Placement(transformation(
          extent={{-10,11},{10,-11}},
          rotation=270,
          origin={-62,-169})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder1(
      r={0,0,0.105*0.5},
      r_shape={0,0,-0.105*0.5},
      diameter=0.16,
      r_0(start=start))
                     annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-44,-238})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation1(r={
          0.216,-0.042,0}) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-74,-52})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation2(r={
          0.224,0.01961,0}) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-74,-18})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation3(r={-0.033266,
          0.156504,0})
      annotation (Placement(transformation(extent={{-94,-54},{-114,-34}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation4(r={
          0.206,0.029,0}) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-96,66})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation5(r={
          0.104,-0.094,0}) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-66,66})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute1(animation=false)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-94,90})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute3(animation=false)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-64,90})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder(
      r={0,0,0.15*0.5},
      r_shape={0,0,-0.15*0.5},
      diameter=0.135) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-98,122})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder2(
      r={0,0,0.15*0.5},
      r_shape={0,0,-0.15*0.5},
      diameter=0.125) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-56,122})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute4(animation=false)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={-136,-44})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder3(
      r={0,0,0.15*0.5},
      r_shape={0,0,-0.15*0.5},
      diameter=0.08)
      annotation (Placement(transformation(extent={{-156,-54},{-176,-34}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation6(r={
          0.17067,0.088845,0})
      annotation (Placement(transformation(extent={{-102,-28},{-122,-8}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation7(r={
          0.065,-0.076,0})
      annotation (Placement(transformation(extent={{-62,-76},{-42,-56}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute5(useAxisFlange=false,
        animation=false)
      annotation (Placement(transformation(extent={{-48,-50},{-28,-30}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute6(useAxisFlange=false,
        animation=false)
      annotation (Placement(transformation(extent={{-48,-22},{-28,-2}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation8(r={
          0.076368,0.076368,0})
      annotation (Placement(transformation(extent={{-10,-40},{10,-20}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation9(r={
          0.115,-0.090,0})
      annotation (Placement(transformation(extent={{-14,-70},{6,-50}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation10(r={-0.076368,
          0.076368,0})
      annotation (Placement(transformation(extent={{-4,16},{16,36}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation11(r={-0.115,
          -0.090,0})
      annotation (Placement(transformation(extent={{-8,-14},{12,6}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute7(useAxisFlange=false,
        animation=false)
      annotation (Placement(transformation(extent={{24,-14},{44,6}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute8(useAxisFlange=false,
        animation=false)
      annotation (Placement(transformation(extent={{18,-70},{38,-50}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation12(r={
          0.051357,-0.016687,0})
      annotation (Placement(transformation(extent={{72,4},{92,24}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation13(r={-0.051357,
          -0.016687,0})
      annotation (Placement(transformation(extent={{72,-20},{92,0}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation14(r={
          0.051357,-0.016687,0})
      annotation (Placement(transformation(extent={{60,-46},{80,-26}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation15(r={-0.051357,
          -0.016687,0})
      annotation (Placement(transformation(extent={{62,-98},{82,-78}})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder4(
      r={0,0,0.023*0.5},
      r_shape={0,0,-0.023*0.5},
      diameter=0.07)
      annotation (Placement(transformation(extent={{108,4},{128,24}})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder5(
      r={0,0,0.023*0.5},
      r_shape={0,0,-0.023*0.5},
      diameter=0.07)
      annotation (Placement(transformation(extent={{108,-18},{128,2}})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder6(
      r={0,0,0.023*0.5},
      r_shape={0,0,-0.023*0.5},
      diameter=0.07)
      annotation (Placement(transformation(extent={{98,-44},{118,-24}})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder7(
      r={0,0,0.023*0.5},
      r_shape={0,0,-0.023*0.5},
      diameter=0.07)
      annotation (Placement(transformation(extent={{106,-82},{126,-62}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation16(r={
          0.0949,0.0373831,0})
      annotation (Placement(transformation(extent={{-14,-138},{6,-118}})));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel
      springDamperParallel(
      c=1000000,
      s_unstretched=0.12,
      d=10000) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={4,-102})));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel
      springDamperParallel1(
      c=1000000,
      s_unstretched=0.15,
      d=10000) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={48,60})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation17(r={
          0.206,0.029,0} + {-0.12,0.04,0}) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-120,26})));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel
      springDamperParallel2(
      c=1000000,
      s_unstretched=0.14,
      d=10000) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-156,6})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute9(useAxisFlange=false,
        animation=false) annotation (Placement(transformation(
          extent={{-7,-7},{7,7}},
          rotation=90,
          origin={-87,3})));
    Modelica.Mechanics.MultiBody.Joints.Prismatic prismatic(useAxisFlange=true)
                       annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-86,30})));
    Modelica.Mechanics.Translational.Sources.Position position(exact=true)
      annotation (Placement(transformation(extent={{-148,56},{-128,76}})));
    Modelica.Blocks.Sources.RealExpression realExpression
      annotation (Placement(transformation(extent={{-192,56},{-172,76}})));
    contact contact1(R=0.07*0.5)
      annotation (Placement(transformation(extent={{156,20},{176,40}})));
    contact contact2(R=0.07*0.5)
      annotation (Placement(transformation(extent={{156,-18},{176,2}})));
    contact contact3(R=0.07*0.5)
      annotation (Placement(transformation(extent={{152,-50},{172,-30}})));
    contact contact4(R=0.07*0.5)
      annotation (Placement(transformation(extent={{144,-92},{164,-72}})));
    contact contact5(R=0.125*0.5)
      annotation (Placement(transformation(extent={{-12,96},{8,116}})));
    contact contact6(R=0.135*0.5)
      annotation (Placement(transformation(extent={{-44,152},{-24,172}})));
    contact contact7(R=0.08*0.5)
      annotation (Placement(transformation(extent={{-192,-40},{-212,-20}})));
    contact_sprocket contact_sprocket1
      annotation (Placement(transformation(extent={{102,-196},{122,-176}})));
    rotational_stop rotational_stop1
      annotation (Placement(transformation(extent={{30,12},{46,26}})));
    rotational_stop rotational_stop2
      annotation (Placement(transformation(extent={{24,-44},{40,-30}})));
     Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a annotation (
      Placement(visible = true, transformation(origin={-256,-68},    extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {-100, -2}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b[80] annotation (
      Placement(visible = true, transformation(origin = {494, 64}, extent = {{-16, -16}, {16, 16}}, rotation = 0), iconTransformation(origin = {98, 0}, extent = {{-16, -16}, {16, 16}}, rotation = 0)));
      Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_a annotation (
        Placement(transformation(extent={{-10,202},{10,222}}),
                                                             iconTransformation(
            extent={{-10,202},{10,222}})));
    track_with_grouser cinghiafinal(
      start=start,
      dh=dh,
      l=0.0586*38,
      startX={0.23480135,0.26272747,0.29065832,0.3185929,0.34652826,0.3744644,
          0.40240145,0.43033797,0.4582741,0.48620963,0.51414156,0.54203176,
          0.5691877,0.592494,0.6125258,0.6317475,0.6507775,0.6697203,0.68844366,
          0.7061789,0.7188599,0.72164524,0.71400744,0.69707525,0.67302054,
          0.6460454,0.6186616,0.5912056,0.5637359,0.5362644,0.5087938,
          0.48132443,0.45385656,0.42639017,0.39892524,0.37146184,0.34399983,
          0.3165391,0.28907877,0.26161474,0.23412834,0.20654036,0.17864704,
          0.15156525,0.12551455,0.09973547,0.07401498,0.048304882,0.022594681,-0.0031178982,
          -0.028833333,-0.054551702,-0.08027288,-0.105996236,-0.13171875,-0.15742654,
          -0.1830544,-0.20828277,-0.23124328,-0.24828674,-0.25737712,-0.2576955,
          -0.24921042,-0.23278509,-0.21002477,-0.18335265,-0.1556219,-0.12777959,
          -0.099918306,-0.07205272,-0.044185717,-0.016317738,0.011551864,
          0.03942658,0.06732078,0.09522102,0.12312411,0.15103368,0.17895454,
          0.20687765} + fill(0.1749185, 80),
      startY={-0.19050065,-0.19128436,-0.191876,-0.19225277,-0.19256467,-0.19280101,
          -0.19275837,-0.19257312,-0.19234018,-0.19204272,-0.19150436,-0.18988861,
          -0.18333949,-0.16794021,-0.14846705,-0.12819362,-0.1077401,-0.08720573,
          -0.066471115,-0.04488537,-0.019996235,0.0077972948,0.034665816,
          0.056881763,0.07108041,0.07834999,0.08388642,0.08905316,0.094146565,
          0.09923051,0.10431874,0.1094142,0.11451752,0.11962885,0.12474819,
          0.12987551,0.1350106,0.14015242,0.14529632,0.15042046,0.15542334,
          0.15983148,0.16139533,0.15453944,0.1444472,0.13367996,0.12277352,
          0.11184263,0.100912035,0.089987054,0.079068825,0.068157546,
          0.057252925,0.046353456,0.03545206,0.024515986,0.013393934,
          0.0013935236,-0.014518096,-0.03665005,-0.06306362,-0.09099584,-0.11760997,
          -0.14020458,-0.15639916,-0.16470166,-0.1680903,-0.17039005,-0.17244737,
          -0.1744456,-0.17642388,-0.17838836,-0.18032962,-0.18219599,-0.18374395,
          -0.18517882,-0.18655726,-0.18779744,-0.18875031,-0.18963498} + fill(
          0.0689022, 80),
      startZ=fill(0, 80),
      thZ(displayUnit="rad") = {-0.028048225,-0.021213358,-0.0134947635,-0.011164346,
        -0.008451123,0.0015039145,0.0066258255,0.00833827,0.010655379,
        0.019309895,0.058046427,0.23605531,0.58301604,0.7710687,0.8119807,
        0.82142866,0.8257008,0.8363717,0.88323003,1.0989887,1.4702258,1.8471031,
        2.221385,2.6072214,2.878054,2.9420385,2.9555721,2.9582531,2.9586,
        2.9584448,2.9581816,2.9578953,2.957604,2.957312,2.9570217,2.9567392,
        2.9564936,2.9564195,2.9571424,2.9615734,2.9832513,3.0851452,-2.8942156,
        -2.772112,-2.74597,-2.7405431,-2.7395878,-2.7395978,-2.739816,-2.7400784,
        -2.740348,-2.740607,-2.7408068,-2.7407305,-2.7393756,-2.7321024,-2.6974335,
        -2.5360084,-2.2276287,-1.9028283,-1.5827664,-1.2627372,-0.94276667,-0.61906016,
        -0.30260813,-0.121777736,-0.08245181,-0.073716834,-0.07158942,-0.070871994,
        -0.0703758,-0.069541104,-0.066845916,-0.05545373,-0.0513865,-0.049357902,
        -0.044453796,-0.034124646,-0.03167472,-0.030992724},
      there_is_grouser={true,true,true,true,true,true,true,false,true,false,
          true,false,true,true,true,true,true,true,true,true,true,true,true,
          true,true,true,true,true,false,true,false,true,false,true,true,true,
          true,true,true,true,true,true,true,true,true,true,true,false,true,
          false,true,false,true,true,true,true,true,true,true,true,true,true,
          true,true,true,true,true,true,false,true,false,true,false,true,true,
          true,true,true,true,true},
      cm={0.0052,0.0241,0.006665,0.025565,0.00813,0.02703,0.009595,0,0.01172,0,
          0.01465,0,0.0175075,7.25e-05,0.0189725,0.0015375,0.0204375,0.0030025,
          0.0219025,0.0044675,0.0233675,0.0059325,0.0248325,0.0073975,0.0262975,
          0.0088625,0.0277625,0.0103275,0,0.013185,0,0.016115,0,0.01824,
          0.000805,0.019705,0.00227,0.02117,0.003735,0.022635,0.0052,0.0241,
          0.006665,0.025565,0.00813,0.02703,0.009595,0,0.01172,0,0.01465,0,
          0.0175075,7.25e-05,0.0189725,0.0015375,0.0204375,0.0030025,0.0219025,
          0.0044675,0.0233675,0.0059325,0.0248325,0.0073975,0.0262975,0.0088625,
          0.0277625,0.0103275,0,0.013185,0,0.016115,0,0.01824,0.000805,0.019705,
          0.00227,0.02117,0.003735,0.022635},
      lg={0.0104,0.00747,0.01333,0.00454,0.01626,0.00161,0.01919,0,0.0208,0,
          0.0208,0,0.020655,0.000145,0.017725,0.003075,0.014795,0.006005,
          0.011865,0.008935,0.008935,0.011865,0.006005,0.014795,0.003075,
          0.017725,0.000145,0.020655,0,0.0208,0,0.0208,0,0.01919,0.00161,
          0.01626,0.00454,0.01333,0.00747,0.0104,0.0104,0.00747,0.01333,0.00454,
          0.01626,0.00161,0.01919,0,0.0208,0,0.0208,0,0.020655,0.000145,
          0.017725,0.003075,0.014795,0.006005,0.011865,0.008935,0.008935,
          0.011865,0.006005,0.014795,0.003075,0.017725,0.000145,0.020655,0,
          0.0208,0,0.0208,0,0.01919,0.00161,0.01626,0.00454,0.01333,0.00747,
          0.0104}) annotation (Placement(visible=true, transformation(
          origin={315,11},
          extent={{-67,-67},{67,67}},
          rotation=0)));
    Modelica.Mechanics.Rotational.Components.Damper damper(d=7)
      annotation (Placement(transformation(extent={{10,10},{-10,-10}},
          rotation=90,
          origin={-94,-166})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation18(r={
          0.3,0,0})
      annotation (Placement(transformation(extent={{-136,-92},{-156,-72}})));
  equation
    connect(revolute.frame_b, fixedTranslation.frame_a) annotation (Line(
        points={{-44,-80},{-38,-80},{-38,-124},{-50,-124}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute2.frame_a, fixedTranslation.frame_b) annotation (Line(
        points={{-62,-159},{-62,-150},{-50,-150},{-50,-144}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute2.frame_b, bodyCylinder1.frame_a) annotation (Line(
        points={{-62,-179},{-62,-208},{-44,-208},{-44,-228}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation2.frame_a, fixedTranslation1.frame_b) annotation (
        Line(
        points={{-74,-28},{-74,-42}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation3.frame_a, fixedTranslation1.frame_b) annotation (
        Line(
        points={{-94,-44},{-88,-44},{-88,-34},{-74,-34},{-74,-42}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute1.frame_a, fixedTranslation4.frame_b) annotation (Line(
        points={{-94,80},{-94,76},{-96,76}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute3.frame_a, fixedTranslation5.frame_b) annotation (Line(
        points={{-64,80},{-64,76},{-66,76}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder.frame_a, revolute1.frame_b) annotation (Line(
        points={{-98,112},{-98,106},{-94,106},{-94,100}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder2.frame_a, revolute3.frame_b) annotation (Line(
        points={{-56,112},{-56,106},{-64,106},{-64,100}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute4.frame_a, fixedTranslation3.frame_b) annotation (Line(
        points={{-126,-44},{-114,-44}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder3.frame_a, revolute4.frame_b) annotation (Line(
        points={{-156,-44},{-146,-44}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation6.frame_a, fixedTranslation1.frame_b) annotation (
        Line(
        points={{-102,-18},{-90,-18},{-90,-30},{-74,-30},{-74,-42}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute6.frame_a, fixedTranslation1.frame_b) annotation (Line(
        points={{-48,-12},{-58,-12},{-58,-34},{-74,-34},{-74,-42}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute5.frame_a, fixedTranslation1.frame_b) annotation (Line(
        points={{-48,-40},{-60,-40},{-60,-34},{-74,-34},{-74,-42}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation8.frame_a, revolute5.frame_b) annotation (Line(
        points={{-10,-30},{-18,-30},{-18,-38},{-28,-38},{-28,-40}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute5.frame_b, fixedTranslation9.frame_a) annotation (Line(
        points={{-28,-40},{-28,-60},{-14,-60}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation11.frame_a, revolute6.frame_b) annotation (Line(
        points={{-8,-4},{-22,-4},{-22,-12},{-28,-12}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation10.frame_a, revolute6.frame_b) annotation (Line(
        points={{-4,26},{-16,26},{-16,22},{-28,22},{-28,-12}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation11.frame_b, revolute7.frame_a) annotation (Line(
        points={{12,-4},{24,-4}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation9.frame_b, revolute8.frame_a) annotation (Line(
        points={{6,-60},{18,-60}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation12.frame_a, revolute7.frame_b) annotation (Line(
        points={{72,14},{50,14},{50,-4},{44,-4}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute7.frame_b, fixedTranslation13.frame_a) annotation (Line(
        points={{44,-4},{66,-4},{66,-10},{72,-10}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation14.frame_a, revolute8.frame_b) annotation (Line(
        points={{60,-36},{44,-36},{44,-60},{38,-60}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute8.frame_b, fixedTranslation15.frame_a) annotation (Line(
        points={{38,-60},{38,-80},{62,-80},{62,-88}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation12.frame_b, bodyCylinder4.frame_a) annotation (Line(
        points={{92,14},{108,14}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation13.frame_b, bodyCylinder5.frame_a) annotation (Line(
        points={{92,-10},{92,-8},{108,-8}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation14.frame_b, bodyCylinder6.frame_a) annotation (Line(
        points={{80,-36},{80,-34},{98,-34}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder7.frame_a, fixedTranslation15.frame_b) annotation (Line(
        points={{106,-72},{86,-72},{86,-88},{82,-88}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation.frame_a, fixedTranslation16.frame_a) annotation (
        Line(
        points={{-50,-124},{-50,-120},{-14,-120},{-14,-128}},
        color={95,95,95},
        thickness=0.5));
    connect(springDamperParallel.frame_a, fixedTranslation16.frame_b)
      annotation (Line(
        points={{4,-112},{6,-112},{6,-128}},
        color={95,95,95},
        thickness=0.5));
    connect(springDamperParallel.frame_b, fixedTranslation7.frame_b)
      annotation (Line(
        points={{4,-92},{4,-74},{-36,-74},{-36,-66},{-42,-66}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation8.frame_b, springDamperParallel1.frame_a)
      annotation (Line(
        points={{10,-30},{26,-30},{26,36},{48,36},{48,50}},
        color={95,95,95},
        thickness=0.5));
    connect(springDamperParallel1.frame_b, fixedTranslation10.frame_b)
      annotation (Line(
        points={{48,70},{36,70},{36,66},{16,66},{16,26}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation6.frame_b, springDamperParallel2.frame_a)
      annotation (Line(
        points={{-122,-18},{-138,-18},{-138,-4},{-156,-4}},
        color={95,95,95},
        thickness=0.5));
    connect(springDamperParallel2.frame_b, fixedTranslation17.frame_b)
      annotation (Line(
        points={{-156,16},{-152,16},{-152,26},{-120,26},{-120,16}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute9.frame_a, fixedTranslation2.frame_b) annotation (Line(
        points={{-87,-4},{-82,-4},{-82,-8},{-74,-8}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation17.frame_a, revolute9.frame_b) annotation (Line(
        points={{-120,36},{-112,36},{-112,10},{-87,10}},
        color={95,95,95},
        thickness=0.5));
    connect(prismatic.frame_a, revolute9.frame_b) annotation (Line(
        points={{-86,20},{-86,15},{-87,15},{-87,10}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation4.frame_a, prismatic.frame_b) annotation (Line(
        points={{-96,56},{-96,40},{-86,40}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation5.frame_a, prismatic.frame_b) annotation (Line(
        points={{-66,56},{-66,46},{-86,46},{-86,40}},
        color={95,95,95},
        thickness=0.5));
    connect(position.flange, prismatic.axis) annotation (Line(points={{-128,66},
            {-124,66},{-124,50},{-92,50},{-92,38}}, color={0,127,0}));
    connect(realExpression.y, position.s_ref)
      annotation (Line(points={{-171,66},{-150,66}}, color={0,0,127}));
    connect(contact1.frame_b, fixedTranslation12.frame_b) annotation (Line(
        points={{154.2,30.4},{92,30.4},{92,14}},
        color={95,95,95},
        thickness=0.5));
    connect(contact2.frame_b, fixedTranslation13.frame_b) annotation (Line(
        points={{154.2,-7.6},{144,-7.6},{144,-22},{92,-22},{92,-10}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation14.frame_b, contact3.frame_b) annotation (Line(
        points={{80,-36},{82,-36},{82,-58},{150.2,-58},{150.2,-39.6}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation15.frame_b, contact4.frame_b) annotation (Line(
        points={{82,-88},{100,-88},{100,-100},{142.2,-100},{142.2,-81.6}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute3.frame_b, contact5.frame_b) annotation (Line(
        points={{-64,100},{-46,100},{-46,104},{-13.8,104},{-13.8,106.4}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute1.frame_b, contact6.frame_b) annotation (Line(
        points={{-94,100},{-94,162.4},{-45.8,162.4}},
        color={95,95,95},
        thickness=0.5));
    connect(contact7.frame_b, revolute4.frame_b) annotation (Line(
        points={{-190.2,-29.6},{-146,-29.6},{-146,-44}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute2.frame_b, contact_sprocket1.frame_b) annotation (Line(
        points={{-62,-179},{-2,-179},{-2,-216},{100.2,-216},{100.2,-185.6}},
        color={95,95,95},
        thickness=0.5));
    connect(rotational_stop1.frame_a, revolute7.frame_b) annotation (Line(
        points={{46,19},{50,19},{50,-4},{44,-4}},
        color={95,95,95},
        thickness=0.5));
    connect(rotational_stop1.frame_b, revolute7.frame_a) annotation (Line(
        points={{30,19},{24,19},{24,-4}},
        color={95,95,95},
        thickness=0.5));
    connect(rotational_stop2.frame_a, revolute8.frame_b) annotation (Line(
        points={{40,-37},{42,-37},{42,-60},{38,-60}},
        color={95,95,95},
        thickness=0.5));
    connect(rotational_stop2.frame_b, revolute8.frame_a) annotation (Line(
        points={{24,-37},{20,-37},{20,-38},{18,-38},{18,-60}},
        color={95,95,95},
        thickness=0.5));

    connect(contact1.frame_a, cinghiafinal.frame_b) annotation (Line(
        points={{176,30.2},{206,30.2},{206,36},{248,36},{248,9.66}},
        color={95,95,95},
        thickness=0.5));
    connect(contact2.frame_a, cinghiafinal.frame_b) annotation (Line(
        points={{176,-7.8},{212,-7.8},{212,9.66},{248,9.66}},
        color={95,95,95},
        thickness=0.5));
    connect(contact3.frame_a, cinghiafinal.frame_b) annotation (Line(
        points={{172,-39.8},{210,-39.8},{210,-24},{248,-24},{248,9.66}},
        color={95,95,95},
        thickness=0.5));
    connect(contact4.frame_a, cinghiafinal.frame_b) annotation (Line(
        points={{164,-81.8},{194,-81.8},{194,-64},{248,-64},{248,9.66}},
        color={95,95,95},
        thickness=0.5));
    connect(contact_sprocket1.frame_a, cinghiafinal.frame_b) annotation (Line(
        points={{122,-185.8},{162,-185.8},{162,-164},{248,-164},{248,9.66}},
        color={95,95,95},
        thickness=0.5));
    connect(contact5.frame_a, cinghiafinal.frame_b) annotation (Line(
        points={{8,106.2},{74,106.2},{74,102},{248,102},{248,9.66}},
        color={95,95,95},
        thickness=0.5));
    connect(contact6.frame_a, cinghiafinal.frame_b) annotation (Line(
        points={{-24,162.2},{106,162.2},{106,138},{248,138},{248,9.66}},
        color={95,95,95},
        thickness=0.5));
    connect(contact7.frame_a, cinghiafinal.frame_b) annotation (Line(
        points={{-212,-29.8},{-212,176},{248,176},{248,9.66}},
        color={95,95,95},
        thickness=0.5));
    connect(cinghiafinal.frame_a, frame_b) annotation (Line(
        points={{382,5.64},{418,5.64},{418,18},{494,18},{494,64}},
        color={95,95,95},
        thickness=0.5));
    connect(flange_a, revolute2.axis) annotation (Line(points={{0,212},{-112,212},
            {-112,198},{-304,198},{-304,-169},{-73,-169}}, color={0,0,0}));
    connect(damper.flange_a, revolute2.support) annotation (Line(points={{-94,
            -156},{-73,-156},{-73,-163}}, color={0,0,0}));
    connect(damper.flange_b, revolute2.axis) annotation (Line(points={{-94,-176},
            {-86,-176},{-86,-190},{-73,-190},{-73,-169}}, color={0,0,0}));

    connect(fixedTranslation1.frame_a, fixedTranslation18.frame_a) annotation (
        Line(
        points={{-74,-62},{-78,-62},{-78,-84},{-136,-84},{-136,-82}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation7.frame_a, fixedTranslation18.frame_a) annotation (
        Line(
        points={{-62,-66},{-130,-66},{-130,-70},{-136,-70},{-136,-82}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute.frame_a, fixedTranslation18.frame_a) annotation (Line(
        points={{-64,-80},{-84,-80},{-84,-106},{-136,-106},{-136,-82}},
        color={95,95,95},
        thickness=0.5));
    connect(frame_a, fixedTranslation18.frame_b) annotation (Line(
        points={{-256,-68},{-226,-68},{-226,-84},{-156,-84},{-156,-82}},
        color={95,95,95},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end trackterrainfinal1_grouser_pro;

  model chain
    Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a annotation (Placement(
          transformation(extent={{-118,-16},{-86,16}}), iconTransformation(extent=
             {{-118,-16},{-86,16}})));
    Modelica.Mechanics.MultiBody.Interfaces.Frame_b frame_b annotation (Placement(
          transformation(extent={{82,-16},{114,16}}), iconTransformation(extent={{
              82,-16},{114,16}})));
    Modelica.Units.SI.Angle rth[3];
    Modelica.Mechanics.MultiBody.Frames.Orientation Rrel;
  equation
    frame_a.f = {0, 0, 0};
    frame_a.f+frame_b.f = {0, 0, 0};
    frame_a.t = -[0,0,0; 0,0,0;0,0,1e5]*rth - [0,0,0; 0,0,0;0,0,5]*Rrel.w;
    frame_a.t + frame_b.t = {0, 0, 0};
    Rrel = Modelica.Mechanics.MultiBody.Frames.relativeRotation(frame_a.R, frame_b.R);
    der(rth) = Rrel.w;
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end chain;

  model track_bunker_pro_new
    extends Modelica.Icons.Example;
    Modelica.Mechanics.MultiBody.Parts.Fixed fixed
      annotation (Placement(transformation(extent={{-124,-58},{-104,-38}})));
    cinghia cinghia1(
      sX=1.130312455 - 0.2,
      sZ=0,
      l=0.0586*38)
      annotation (Placement(transformation(extent={{380,-24},{432,28}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute(useAxisFlange=true,
        animation=false)
      annotation (Placement(transformation(extent={{-38,-54},{-18,-34}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation(r={-0.18554,
          -0.07342,0})   annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-30,-82})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder1(
      r={0,0,0},
      r_shape={0,0,-0.105*0.5},
      length=0.105,
      diameter=0.16,
      r_0(fixed=false),
      useQuaternions=false)
                     annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-22,-192})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation1(r={
          0.21488,-0.04227,0})
                           annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-52,-6})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation2(r={
          0.22252,0.01926,0})
                            annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-52,28})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation3(r={-0.03379,
          0.15530,0})
      annotation (Placement(transformation(extent={{-72,-8},{-92,12}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation4(r={
          0.20314,0.02838,0})
                          annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-74,112})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation5(r={
          0.10492,-0.09153,0})
                           annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-44,112})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute1(animation=false)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-72,136})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute3(animation=false)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-42,136})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder(
      r={0,0,0},
      r_shape={0,0,-0.15*0.5},
      length=0.15,
      diameter=0.139,
      r_0(fixed=false),
      useQuaternions=false)
                      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-76,168})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder2(
      r={0,0,0},
      r_shape={0,0,-0.15*0.5},
      length=0.15,
      diameter=0.126,
      useQuaternions=false)
                      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-34,168})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute4(animation=false)
      annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=180,
          origin={-114,2})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder3(
      r={0,0,0},
      r_shape={0,0,-0.15*0.5},
      length=0.15,
      diameter=0.08,
      useQuaternions=false)
      annotation (Placement(transformation(extent={{-134,-8},{-154,12}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation6(r={
          0.17068,0.08769,0})
      annotation (Placement(transformation(extent={{-80,18},{-100,38}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation7(r={
          0.06426,-0.07476,0})
      annotation (Placement(transformation(extent={{-40,-32},{-20,-12}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute5(useAxisFlange=true,
        animation=false)
      annotation (Placement(transformation(extent={{-26,-4},{-6,16}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute6(useAxisFlange=true,
        animation=false)
      annotation (Placement(transformation(extent={{-26,24},{-6,44}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation8(r={
          0.07501,0.07623,0})
      annotation (Placement(transformation(extent={{12,6},{32,26}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation9(r={
          0.11355,-0.09018,0})
      annotation (Placement(transformation(extent={{8,-24},{28,-4}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation10(r={-0.07501,
          0.07623,0})
      annotation (Placement(transformation(extent={{18,62},{38,82}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation11(r={-0.11355,
          -0.09018,0})
      annotation (Placement(transformation(extent={{14,32},{34,52}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute7(useAxisFlange=true,
        animation=false)
      annotation (Placement(transformation(extent={{46,32},{66,52}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute8(useAxisFlange=true,
        animation=false)
      annotation (Placement(transformation(extent={{40,-24},{60,-4}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation12(r={
          0.05250,-0.01633,0})
      annotation (Placement(transformation(extent={{94,50},{114,70}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation13(r={-0.05250,
          -0.01633,0})
      annotation (Placement(transformation(extent={{94,26},{114,46}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation14(r={
          0.05250,-0.01633,0})
      annotation (Placement(transformation(extent={{82,0},{102,20}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation15(r={-0.05250,
          -0.01633,0})
      annotation (Placement(transformation(extent={{84,-52},{104,-32}})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder4(
      r={0,0,0},
      r_shape={0,0,-0.023*0.5},
      length=0.023,
      diameter=0.069,
      angles_fixed=false,
      w_0_fixed=false,
      useQuaternions=false)
      annotation (Placement(transformation(extent={{156,48},{176,68}})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder5(
      r={0,0,0},
      r_shape={0,0,-0.023*0.5},
      length=0.023,
      diameter=0.069,
      angles_fixed=false,
      w_0_fixed=false,
      useQuaternions=false)
      annotation (Placement(transformation(extent={{150,30},{170,50}})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder6(
      r={0,0,0},
      r_shape={0,0,-0.023*0.5},
      length=0.023,
      diameter=0.069,
      angles_fixed=false,
      w_0_fixed=false,
      useQuaternions=false)
      annotation (Placement(transformation(extent={{144,0},{164,20}})));
    Modelica.Mechanics.MultiBody.Parts.BodyCylinder bodyCylinder7(
      r={0,0,0},
      r_shape={0,0,-0.023*0.5},
      length=0.023,
      diameter=0.069,
      angles_fixed=false,
      w_0_fixed=false,
      useQuaternions=false)
      annotation (Placement(transformation(extent={{140,-36},{160,-16}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation16(r={
          0.084035,0.0331024,0})
      annotation (Placement(transformation(extent={{8,-92},{28,-72}})));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel
      springDamperParallel(
      c=100000,
      s_unstretched=0.10,
      d=10000)
              annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={26,-56})));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel
      springDamperParallel1(
      c=100000,
      s_unstretched=0.16,
      d=10000)
              annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={70,106})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation17(r={
          0.20314,0.02838,0} + {-0.12,0.04005,0})
                                           annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-98,72})));
    Modelica.Mechanics.MultiBody.Forces.SpringDamperParallel
      springDamperParallel2(
      c=100000,
      s_unstretched=0.13,
      d=10000)
              annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-134,52})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute9(useAxisFlange=true,
        animation=false) annotation (Placement(transformation(
          extent={{-7,-7},{7,7}},
          rotation=90,
          origin={-65,49})));
    Modelica.Mechanics.MultiBody.Joints.Prismatic prismatic(useAxisFlange=true,
        animation=false)
                       annotation (Placement(transformation(
          extent={{-10,-10},{10,10}},
          rotation=90,
          origin={-64,76})));
    Modelica.Mechanics.Translational.Sources.Position position(exact=true)
      annotation (Placement(transformation(extent={{-126,102},{-106,122}})));
    Modelica.Blocks.Sources.RealExpression realExpression
      annotation (Placement(transformation(extent={{-170,102},{-150,122}})));
    contact contact1(
      R=0.069*0.5,
      Krx=10,
      Drx=0.1)
      annotation (Placement(transformation(extent={{178,66},{198,86}})));
    contact contact2(
      R=0.069*0.5,
      Krx=10,
      Drx=0.1)
      annotation (Placement(transformation(extent={{178,28},{198,48}})));
    contact contact3(
      R=0.069*0.5,
      Krx=10,
      Drx=0.1)
      annotation (Placement(transformation(extent={{174,-4},{194,16}})));
    contact contact4(
      R=0.069*0.5,
      Krx=10,
      Drx=0.1)
      annotation (Placement(transformation(extent={{166,-46},{186,-26}})));
    contact contact5(R=0.126*0.5)
      annotation (Placement(transformation(extent={{10,142},{30,162}})));
    contact contact6(R=0.139*0.5)
      annotation (Placement(transformation(extent={{-22,198},{-2,218}})));
    contact contact7(R=0.08*0.5)
      annotation (Placement(transformation(extent={{-170,6},{-190,26}})));
    contact_sprocket contact_sprocket1
      annotation (Placement(transformation(extent={{124,-150},{144,-130}})));
    rotational_stop rotational_stop1
      annotation (Placement(transformation(extent={{52,58},{68,72}})));
    rotational_stop rotational_stop2
      annotation (Placement(transformation(extent={{46,2},{62,16}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute10(useAxisFlange=false,
        animation=false)
      annotation (Placement(transformation(extent={{-10,-10},{10,10}},
          rotation=270,
          origin={-12,-138})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute11(useAxisFlange=false,
        animation=false)
      annotation (Placement(transformation(extent={{124,50},{144,70}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute12(useAxisFlange=false,
        animation=false)
      annotation (Placement(transformation(extent={{120,26},{140,46}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute13(useAxisFlange=false,
        animation=false)
      annotation (Placement(transformation(extent={{112,-2},{132,18}})));
    Modelica.Mechanics.MultiBody.Joints.Revolute revolute14(useAxisFlange=false,
        animation=false)
      annotation (Placement(transformation(extent={{110,-36},{130,-16}})));
    Modelica.Mechanics.Rotational.Sources.Position position1(exact=true)
      annotation (Placement(transformation(extent={{-96,-104},{-76,-84}})));
    Modelica.Blocks.Sources.RealExpression realExpression1
      annotation (Placement(transformation(extent={{-150,-104},{-130,-84}})));
    Modelica.Mechanics.Rotational.Sources.Position position2(exact=true)
      annotation (Placement(transformation(extent={{-244,-58},{-224,-38}})));
    Modelica.Blocks.Sources.RealExpression realExpression2
      annotation (Placement(transformation(extent={{-298,-58},{-278,-38}})));
    Modelica.Mechanics.Rotational.Sources.Position position3(exact=true)
      annotation (Placement(transformation(extent={{-268,2},{-248,22}})));
    Modelica.Blocks.Sources.RealExpression realExpression3
      annotation (Placement(transformation(extent={{-322,2},{-302,22}})));
    Modelica.Mechanics.Rotational.Sources.Position position4(exact=true)
      annotation (Placement(transformation(extent={{-226,50},{-206,70}})));
    Modelica.Blocks.Sources.RealExpression realExpression4
      annotation (Placement(transformation(extent={{-280,50},{-260,70}})));
    Modelica.Mechanics.Rotational.Sources.Position position5(exact=true)
      annotation (Placement(transformation(extent={{46,-112},{66,-92}})));
    Modelica.Blocks.Sources.RealExpression realExpression5
      annotation (Placement(transformation(extent={{-8,-112},{12,-92}})));
    Modelica.Mechanics.Rotational.Sources.Position position6(exact=true)
      annotation (Placement(transformation(extent={{134,116},{154,136}})));
    Modelica.Blocks.Sources.RealExpression realExpression6
      annotation (Placement(transformation(extent={{80,116},{100,136}})));
  equation
    connect(revolute.frame_b,fixedTranslation. frame_a) annotation (Line(
        points={{-18,-44},{-12,-44},{-12,-72},{-30,-72}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation2.frame_a,fixedTranslation1. frame_b) annotation (
        Line(
        points={{-52,18},{-52,4}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation3.frame_a,fixedTranslation1. frame_b) annotation (
        Line(
        points={{-72,2},{-66,2},{-66,12},{-52,12},{-52,4}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute1.frame_a,fixedTranslation4. frame_b) annotation (Line(
        points={{-72,126},{-72,122},{-74,122}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute3.frame_a,fixedTranslation5. frame_b) annotation (Line(
        points={{-42,126},{-42,122},{-44,122}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder.frame_a,revolute1. frame_b) annotation (Line(
        points={{-76,158},{-76,152},{-72,152},{-72,146}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder2.frame_a,revolute3. frame_b) annotation (Line(
        points={{-34,158},{-34,152},{-42,152},{-42,146}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute4.frame_a,fixedTranslation3. frame_b) annotation (Line(
        points={{-104,2},{-92,2}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder3.frame_a,revolute4. frame_b) annotation (Line(
        points={{-134,2},{-124,2}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation6.frame_a,fixedTranslation1. frame_b) annotation (
        Line(
        points={{-80,28},{-68,28},{-68,16},{-52,16},{-52,4}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute6.frame_a,fixedTranslation1. frame_b) annotation (Line(
        points={{-26,34},{-36,34},{-36,12},{-52,12},{-52,4}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute5.frame_a,fixedTranslation1. frame_b) annotation (Line(
        points={{-26,6},{-38,6},{-38,12},{-52,12},{-52,4}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation8.frame_a,revolute5. frame_b) annotation (Line(
        points={{12,16},{4,16},{4,8},{-6,8},{-6,6}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute5.frame_b,fixedTranslation9. frame_a) annotation (Line(
        points={{-6,6},{-6,-14},{8,-14}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation11.frame_a,revolute6. frame_b) annotation (Line(
        points={{14,42},{0,42},{0,34},{-6,34}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation10.frame_a,revolute6. frame_b) annotation (Line(
        points={{18,72},{6,72},{6,68},{-6,68},{-6,34}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation11.frame_b,revolute7. frame_a) annotation (Line(
        points={{34,42},{46,42}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation9.frame_b,revolute8. frame_a) annotation (Line(
        points={{28,-14},{40,-14}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation12.frame_a,revolute7. frame_b) annotation (Line(
        points={{94,60},{72,60},{72,42},{66,42}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute7.frame_b,fixedTranslation13. frame_a) annotation (Line(
        points={{66,42},{88,42},{88,36},{94,36}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation14.frame_a,revolute8. frame_b) annotation (Line(
        points={{82,10},{66,10},{66,-14},{60,-14}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute8.frame_b,fixedTranslation15. frame_a) annotation (Line(
        points={{60,-14},{60,-34},{84,-34},{84,-42}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation.frame_a,fixedTranslation16. frame_a) annotation (
        Line(
        points={{-30,-72},{-30,-74},{8,-74},{8,-82}},
        color={95,95,95},
        thickness=0.5));
    connect(springDamperParallel.frame_a,fixedTranslation16. frame_b)
      annotation (Line(
        points={{26,-66},{28,-66},{28,-82}},
        color={95,95,95},
        thickness=0.5));
    connect(springDamperParallel.frame_b,fixedTranslation7. frame_b)
      annotation (Line(
        points={{26,-46},{26,-28},{-14,-28},{-14,-22},{-20,-22}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation8.frame_b,springDamperParallel1. frame_a)
      annotation (Line(
        points={{32,16},{48,16},{48,82},{70,82},{70,96}},
        color={95,95,95},
        thickness=0.5));
    connect(springDamperParallel1.frame_b,fixedTranslation10. frame_b)
      annotation (Line(
        points={{70,116},{58,116},{58,112},{38,112},{38,72}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation6.frame_b,springDamperParallel2. frame_a)
      annotation (Line(
        points={{-100,28},{-116,28},{-116,42},{-134,42}},
        color={95,95,95},
        thickness=0.5));
    connect(springDamperParallel2.frame_b,fixedTranslation17. frame_b)
      annotation (Line(
        points={{-134,62},{-130,62},{-130,72},{-98,72},{-98,62}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute9.frame_a,fixedTranslation2. frame_b) annotation (Line(
        points={{-65,42},{-60,42},{-60,38},{-52,38}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation17.frame_a,revolute9. frame_b) annotation (Line(
        points={{-98,82},{-90,82},{-90,56},{-65,56}},
        color={95,95,95},
        thickness=0.5));
    connect(prismatic.frame_a,revolute9. frame_b) annotation (Line(
        points={{-64,66},{-64,61},{-65,61},{-65,56}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation4.frame_a,prismatic. frame_b) annotation (Line(
        points={{-74,102},{-74,86},{-64,86}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation5.frame_a,prismatic. frame_b) annotation (Line(
        points={{-44,102},{-44,92},{-64,92},{-64,86}},
        color={95,95,95},
        thickness=0.5));
    connect(position.flange,prismatic. axis) annotation (Line(points={{-106,112},
            {-102,112},{-102,96},{-70,96},{-70,84}},color={0,127,0}));
    connect(realExpression.y,position. s_ref)
      annotation (Line(points={{-149,112},{-128,112}},
                                                     color={0,0,127}));
    connect(rotational_stop1.frame_a,revolute7. frame_b) annotation (Line(
        points={{68,65},{72,65},{72,42},{66,42}},
        color={95,95,95},
        thickness=0.5));
    connect(rotational_stop1.frame_b,revolute7. frame_a) annotation (Line(
        points={{52,65},{46,65},{46,42}},
        color={95,95,95},
        thickness=0.5));
    connect(rotational_stop2.frame_a,revolute8. frame_b) annotation (Line(
        points={{62,9},{64,9},{64,-14},{60,-14}},
        color={95,95,95},
        thickness=0.5));
    connect(rotational_stop2.frame_b,revolute8. frame_a) annotation (Line(
        points={{46,9},{42,9},{42,8},{40,8},{40,-14}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation.frame_b,revolute10. frame_a) annotation (Line(
        points={{-30,-92},{-24,-92},{-24,-118},{-12,-118},{-12,-128}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute10.frame_b,bodyCylinder1. frame_a) annotation (Line(
        points={{-12,-148},{-12,-168},{-22,-168},{-22,-182}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation12.frame_b,revolute11. frame_a) annotation (Line(
        points={{114,60},{124,60}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute11.frame_b,bodyCylinder4. frame_a) annotation (Line(
        points={{144,60},{150,60},{150,62},{156,62},{156,58}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation13.frame_b,revolute12. frame_a) annotation (Line(
        points={{114,36},{120,36}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute12.frame_b,bodyCylinder5. frame_a) annotation (Line(
        points={{140,36},{150,36},{150,40}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation14.frame_b,revolute13. frame_a) annotation (Line(
        points={{102,10},{102,8},{112,8}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute13.frame_b,bodyCylinder6. frame_a) annotation (Line(
        points={{132,8},{132,10},{144,10}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation15.frame_b,revolute14. frame_a) annotation (Line(
        points={{104,-42},{106,-42},{106,-26},{110,-26}},
        color={95,95,95},
        thickness=0.5));
    connect(revolute14.frame_b,bodyCylinder7. frame_a) annotation (Line(
        points={{130,-26},{140,-26}},
        color={95,95,95},
        thickness=0.5));
    connect(contact1.frame_b,bodyCylinder4. frame_b) annotation (Line(
        points={{176.2,76.4},{168,76.4},{168,58},{176,58}},
        color={95,95,95},
        thickness=0.5));
    connect(contact2.frame_b,bodyCylinder5. frame_b) annotation (Line(
        points={{176.2,38.4},{176.2,40},{170,40}},
        color={95,95,95},
        thickness=0.5));
    connect(contact3.frame_b,bodyCylinder6. frame_b) annotation (Line(
        points={{172.2,6.4},{164,6.4},{164,10}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder7.frame_b,contact4. frame_b) annotation (Line(
        points={{160,-26},{164.2,-26},{164.2,-35.6}},
        color={95,95,95},
        thickness=0.5));
    connect(contact6.frame_b,bodyCylinder. frame_b) annotation (Line(
        points={{-23.8,208.4},{-76,208.4},{-76,178}},
        color={95,95,95},
        thickness=0.5));
    connect(contact5.frame_b,bodyCylinder2. frame_b) annotation (Line(
        points={{8.2,152.4},{-10,152.4},{-10,186},{-34,186},{-34,178}},
        color={95,95,95},
        thickness=0.5));
    connect(contact7.frame_b,bodyCylinder3. frame_b) annotation (Line(
        points={{-168.2,16.4},{-156,16.4},{-156,2},{-154,2}},
        color={95,95,95},
        thickness=0.5));
    connect(bodyCylinder1.frame_b,contact_sprocket1. frame_b) annotation (Line(
        points={{-22,-202},{26,-202},{26,-200},{88,-200},{88,-139.6},{122.2,
            -139.6}},
        color={95,95,95},
        thickness=0.5));
    connect(fixed.frame_b, revolute.frame_a) annotation (Line(
        points={{-104,-48},{-48,-48},{-48,-44},{-38,-44}},
        color={95,95,95},
        thickness=0.5));
    connect(contact_sprocket1.frame_a, cinghia1.frame_b) annotation (Line(
        points={{144,-139.8},{262,-139.8},{262,1.48},{380,1.48}},
        color={95,95,95},
        thickness=0.5));
    connect(contact4.frame_a, cinghia1.frame_b) annotation (Line(
        points={{186,-35.8},{283,-35.8},{283,1.48},{380,1.48}},
        color={95,95,95},
        thickness=0.5));
    connect(contact3.frame_a, cinghia1.frame_b) annotation (Line(
        points={{194,6.2},{287,6.2},{287,1.48},{380,1.48}},
        color={95,95,95},
        thickness=0.5));
    connect(contact2.frame_a, cinghia1.frame_b) annotation (Line(
        points={{198,38.2},{290,38.2},{290,1.48},{380,1.48}},
        color={95,95,95},
        thickness=0.5));
    connect(contact5.frame_a, cinghia1.frame_b) annotation (Line(
        points={{30,152.2},{206,152.2},{206,1.48},{380,1.48}},
        color={95,95,95},
        thickness=0.5));
    connect(contact6.frame_a, cinghia1.frame_b) annotation (Line(
        points={{-2,208.2},{380,208.2},{380,1.48}},
        color={95,95,95},
        thickness=0.5));
    connect(contact7.frame_a, cinghia1.frame_b) annotation (Line(
        points={{-190,16.2},{-190,-266},{380,-266},{380,1.48}},
        color={95,95,95},
        thickness=0.5));
    connect(realExpression1.y,position1. phi_ref) annotation (Line(points={{-129,
            -94},{-120,-94},{-120,-96},{-98,-96},{-98,-94}},        color={0,0,
            127}));
    connect(position1.flange, revolute.axis) annotation (Line(points={{-76,-94},
            {-46,-94},{-46,-34},{-28,-34}}, color={0,0,0}));
    connect(realExpression2.y,position2. phi_ref) annotation (Line(points={{-277,
            -48},{-268,-48},{-268,-50},{-246,-50},{-246,-48}},      color={0,0,
            127}));
    connect(position2.flange, revolute5.axis) annotation (Line(points={{-224,
            -48},{-152,-48},{-152,-20},{-16,-20},{-16,16}}, color={0,0,0}));
    connect(realExpression3.y,position3. phi_ref) annotation (Line(points={{-301,12},
            {-292,12},{-292,10},{-270,10},{-270,12}},               color={0,0,
            127}));
    connect(position3.flange, revolute6.axis) annotation (Line(points={{-248,12},
            {-230,12},{-230,40},{-16,40},{-16,44}}, color={0,0,0}));
    connect(realExpression4.y,position4. phi_ref) annotation (Line(points={{-259,60},
            {-250,60},{-250,58},{-228,58},{-228,60}},               color={0,0,
            127}));
    connect(position4.flange, revolute9.axis) annotation (Line(points={{-206,60},
            {-154,60},{-154,52},{-72,52},{-72,49}}, color={0,0,0}));
    connect(realExpression5.y,position5. phi_ref) annotation (Line(points={{13,-102},
            {22,-102},{22,-104},{44,-104},{44,-102}},               color={0,0,
            127}));
    connect(position5.flange, revolute8.axis) annotation (Line(points={{66,-102},
            {86,-102},{86,-58},{50,-58},{50,-4}}, color={0,0,0}));
    connect(realExpression6.y,position6. phi_ref) annotation (Line(points={{101,126},
            {110,126},{110,124},{132,124},{132,126}},               color={0,0,
            127}));
    connect(position6.flange, revolute7.axis) annotation (Line(points={{154,126},
            {170,126},{170,90},{84,90},{84,52},{56,52}}, color={0,0,0}));
    connect(contact1.frame_a, cinghia1.frame_b) annotation (Line(
        points={{198,76.2},{242,76.2},{242,62},{380,62},{380,1.48}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation1.frame_a, fixed.frame_b) annotation (Line(
        points={{-52,-16},{-54,-16},{-54,-48},{-104,-48}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation7.frame_a, fixed.frame_b) annotation (Line(
        points={{-40,-22},{-54,-22},{-54,-48},{-104,-48}},
        color={95,95,95},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=30,
        Tolerance=0.001,
        __Dymola_Algorithm="Dassl"));
  end track_bunker_pro_new;

  model IMU
    parameter Boolean anim = true;
    parameter Modelica.Units.SI.Position start[3];
    Modelica.Mechanics.MultiBody.Interfaces.Frame_a frame_a annotation (Placement(
          transformation(extent={{-114,-14},{-82,18}}), iconTransformation(extent=
             {{-114,-14},{-82,18}})));
          Modelica.Units.SI.Velocity v[3];
          Real a[3](each final unit="g");
          Real wx(final unit="deg/s");
          Real wy(final unit="deg/s");
          Real wz(final unit="deg/s");
    Modelica.Mechanics.MultiBody.Parts.PointMass pointMass(
      animation=false,                                     m=1e-6, r_0(start=
            start))
      annotation (Placement(transformation(extent={{24,-10},{44,10}})));
    Modelica.Mechanics.MultiBody.Sensors.CutForce cutForce(animation=false,
        resolveInFrame=Modelica.Mechanics.MultiBody.Types.ResolveInFrameA.frame_a)
      annotation (Placement(transformation(extent={{-48,-8},{-28,12}})));
  equation

    v = Modelica.Mechanics.MultiBody.Frames.resolve2(frame_a.R, der(frame_a.r_0));
    a =((cutForce.force)/pointMass.m)/Modelica.Constants.g_n;

    wx = (frame_a.R.w[1]/Modelica.Constants.pi)*180;
    wy = (frame_a.R.w[2]/Modelica.Constants.pi)*180;
    wz = (frame_a.R.w[3]/Modelica.Constants.pi)*180;
    connect(frame_a, cutForce.frame_a) annotation (Line(
        points={{-98,2},{-48,2}},
        color={95,95,95},
        thickness=0.5));
    connect(pointMass.frame_a, cutForce.frame_b) annotation (Line(
        points={{34,0},{34,2},{-28,2}},
        color={95,95,95},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)));
  end IMU;

  model Fig_15
    parameter Real speed0_ref[1,2001] = -Modelica.Utilities.Streams.readRealMatrix("Fig_15_speeds.mat","speed0", 1,2001);
    parameter Real speed1_ref[1,2001] = -Modelica.Utilities.Streams.readRealMatrix("Fig_15_speeds.mat","speed1", 1,2001);
    parameter Real t[1,2001] = Modelica.Utilities.Streams.readRealMatrix("Fig_15_speeds.mat","t", 1,2001);
    parameter Real table_speed0[2001,2] = transpose([t; speed0_ref]);
    parameter Real table_speed1[2001,2] = transpose([t; speed1_ref]);
    parameter Modelica.Units.SI.Angle a = 0;
    parameter Modelica.Units.SI.Angle b = -0.162526466768425;
    parameter Modelica.Units.SI.Angle g = 0;
    final parameter Modelica.Mechanics.MultiBody.Frames.Orientation R = Modelica.Mechanics.MultiBody.Frames.axesRotations({1,2,3}, {a,b,g}, zeros(3));
    final parameter Modelica.Units.SI.Position r[3] = Modelica.Mechanics.MultiBody.Frames.resolve1(R, {0,0,0.574});
    parameter Modelica.Units.SI.Position start[3] = {5.411049943885385 + 0.012739630961239, 0.45838,2.581958284866688 + 0.002017940438606};
    final parameter Modelica.Units.SI.Position start1[3] = start + Modelica.Mechanics.MultiBody.Frames.resolve1(R, - {0.56507 - 0.42575 - 0.10525,0.29023,-0.574*0.5}-{0,0,0.574});
    final parameter Modelica.Units.SI.Position start2[3] = start + Modelica.Mechanics.MultiBody.Frames.resolve1(R, {-0.18554,-0.07342,0}-{0.21488,-0.04227,0}-{0.56507 - 0.42575 - 0.10525,0.29023,-0.574*0.5});
    extends Modelica.Icons.Example;
    trackterrain_grouser trackterrain_grouser_R(
      start=start2,
      a=a,
      b=b,
      g=g,
      dh={contact_t_R.contact_terrain[i].dh for i in 1:80})
      annotation (Placement(transformation(extent={{74,-8},{138,16}})));
    Modelica.Mechanics.MultiBody.Sensors.AbsoluteAngles Vicon_angles
      annotation (Placement(visible=true, transformation(
          origin={40,-79.5},
          extent={{-14,-13.5},{14,13.5}},
          rotation=0)));
    contact_t_grouser contact_t_R(
      nz=2,
      there_is_grouser={true,true,true,true,true,false,true,false,true,false,
          true,false,true,false,true,false,true,true,true,true,true,true,true,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,true,true,true,true,true,true,true,true,true,false,true,false,
          true,false,true,false,true,false,true,false,true,true,true,true,true,
          true,true,true,false,true,false,true,false,true,false,true,false,true,
          false,true,true,true,true,true},
      cm1={0.0175675,0.0117325,0.0190325,0.0131975,0.0204975,0,0.019573828,0,
          0.01634164,0,0.013109453,0,0.009877265,0,0.0066450783,0,0.00807,
          0.01537,0.009535,0.016835,0.011,0.0183,0.012465,0.019765,0,
          0.021189922,0,0.017957734,0,0.014725547,0,0.01149336,0,0.008261172,0,
          0.0073375,0.0146375,0.0088025,0.0161025,0.0102675,0.0175675,0.0117325,
          0.0190325,0.0131975,0.0204975,0,0.019573828,0,0.01634164,0,
          0.013109453,0,0.009877265,0,0.0066450783,0,0.00807,0.01537,0.009535,
          0.016835,0.011,0.0183,0.012465,0.019765,0,0.021189922,0,0.017957734,0,
          0.014725547,0,0.01149336,0,0.008261172,0,0.0073375,0.0146375,
          0.0088025,0.0161025,0.0102675},
      cm2={0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,0.01172,0,
          0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,0.0029175,
          0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,0.010255,0,
          0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,0.002185,
          0.024185,0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,
          0.01172,0,0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,
          0.0029175,0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,
          0.010255,0,0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,
          0.002185,0.024185},
      lg={0.0073,0.00437,0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,
          0.0146,0,0.0146,0,0.011695,0.002905,0.008765,0.005835,0.005835,
          0.008765,0.002905,0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,
          0.0146,0,0.01316,0.00144,0.01023,0.00437,0.0073,0.0073,0.00437,
          0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,
          0,0.011695,0.002905,0.008765,0.005835,0.005835,0.008765,0.002905,
          0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.01316,
          0.00144,0.01023,0.00437,0.0073},
      fr={true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true,false,true,false,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,true,false,true,false,false,true,false,true,
          false,true,false,true,false,true,false,true,false,true,false,true,
          false,true,false,false,true,false,true,false},
      fl={false,true,false,true,false,false,true,false,true,false,true,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          false,true,false,true,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true,false,false,true,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true})
      annotation (Placement(transformation(extent={{172,-16},{210,22}})));

    contact_t_grouser contact_t_L(
      nz=2,
      there_is_grouser={true,true,true,true,true,false,true,false,true,false,
          true,false,true,false,true,false,true,true,true,true,true,true,true,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,true,true,true,true,true,true,true,true,true,false,true,false,
          true,false,true,false,true,false,true,false,true,true,true,true,true,
          true,true,true,false,true,false,true,false,true,false,true,false,true,
          false,true,true,true,true,true},
      cm1={0.0175675,0.0117325,0.0190325,0.0131975,0.0204975,0,0.019573828,0,
          0.01634164,0,0.013109453,0,0.009877265,0,0.0066450783,0,0.00807,
          0.01537,0.009535,0.016835,0.011,0.0183,0.012465,0.019765,0,
          0.021189922,0,0.017957734,0,0.014725547,0,0.01149336,0,0.008261172,0,
          0.0073375,0.0146375,0.0088025,0.0161025,0.0102675,0.0175675,0.0117325,
          0.0190325,0.0131975,0.0204975,0,0.019573828,0,0.01634164,0,
          0.013109453,0,0.009877265,0,0.0066450783,0,0.00807,0.01537,0.009535,
          0.016835,0.011,0.0183,0.012465,0.019765,0,0.021189922,0,0.017957734,0,
          0.014725547,0,0.01149336,0,0.008261172,0,0.0073375,0.0146375,
          0.0088025,0.0161025,0.0102675},
      cm2={0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,0.01172,0,
          0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,0.0029175,
          0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,0.010255,0,
          0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,0.002185,
          0.024185,0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,
          0.01172,0,0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,
          0.0029175,0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,
          0.010255,0,0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,
          0.002185,0.024185},
      lg={0.0073,0.00437,0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,
          0.0146,0,0.0146,0,0.011695,0.002905,0.008765,0.005835,0.005835,
          0.008765,0.002905,0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,
          0.0146,0,0.01316,0.00144,0.01023,0.00437,0.0073,0.0073,0.00437,
          0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,
          0,0.011695,0.002905,0.008765,0.005835,0.005835,0.008765,0.002905,
          0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.01316,
          0.00144,0.01023,0.00437,0.0073},
      fr={true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true,false,true,false,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,true,false,true,false,false,true,false,true,
          false,true,false,true,false,true,false,true,false,true,false,true,
          false,true,false,false,true,false,true,false},
      fl={false,true,false,true,false,false,true,false,true,false,true,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          false,true,false,true,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true,false,false,true,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true})
      annotation (Placement(transformation(extent={{-128,-14},{-166,24}})));

    inner Modelica.Mechanics.MultiBody.World world
      annotation (Placement(transformation(extent={{-100,-60},{-68,-28}})));
    trackterrain_grouser trackterrain_grouser_L(
      start=start2 - r,
      a=a,
      b=b,
      g=g,
      dh={contact_t_L.contact_terrain[i].dh for i in 1:80})
      annotation (Placement(transformation(extent={{-32,-2},{-96,22}})));
    Modelica.Mechanics.Rotational.Sources.Speed speed_L
      annotation (Placement(transformation(extent={{-108,48},{-80,76}})));
    Modelica.Mechanics.Rotational.Sources.Speed speed_R
      annotation (Placement(transformation(extent={{48,50},{76,78}})));
    Modelica.Mechanics.MultiBody.Parts.BodyShape Chassis(
      animation=false,
      animateSphere=false,
      r={0,0,0.574},
      r_CM={0.03535770079,0.12598,0.574*0.5},
      m=149.4571417177305,
      I_11=(1/12)*149.4571417177305*(0.424^2 + 0.3285^2),
      I_22=(1/12)*149.4571417177305*(0.424^2 + (0.5*0.820)^2) +
          110.347597567632*0.107287700790000^2 + 39.1095441500988*
          0.302712299210000^2,
      I_33=(1/12)*149.4571417177305*((0.5*0.820)^2 + 0.3285^2) +
          110.347597567632*0.107287700790000^2 + 39.1095441500988*
          0.302712299210000^2,
      r_0(fixed=true, start=start1),
      angles_fixed=true,
      angles_start={a,b,g},
      shapeType="box",
      r_shape={0.13307,0.12598,0.075},
      lengthDirection(displayUnit="1") = {0,0,1},
      length=0.424,
      width=0.3285,
      height=0.820,
      useQuaternions=false)
      annotation (Placement(transformation(extent={{-4,-6},{28,26}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation(animation
        =false,                                                          r={
          0.56507 - 0.42575 - 0.10525,0.29023,-0.574*0.5})
      annotation (Placement(transformation(extent={{40,-50},{60,-30}})));
    inner Surface_base Surface(
      nu=2,
      nv=2,
      multiColoredSurface=false,
      color={0,255,0},
      nu1=2,
      nv1=2,
      LengthX=14.8,
      LengthZ=7.3,
      OffsetX=-2.5,
      OffsetZ=-1.7)
                 annotation (Placement(visible=true, transformation(
          origin={147,72.5},
          extent={{-15,-15.5},{15,15.5}},
          rotation=0)));
    Modelica.Blocks.Sources.CombiTimeTable combiTimeTable_L(table=table_speed1,
        timeEvents=Modelica.Blocks.Types.TimeEvents.NoTimeEvents)
      annotation (Placement(transformation(extent={{-170,48},{-140,78}})));
    Modelica.Blocks.Sources.CombiTimeTable combiTimeTable_R(table=table_speed0,
        timeEvents=Modelica.Blocks.Types.TimeEvents.NoTimeEvents)
      annotation (Placement(transformation(extent={{-16,48},{14,78}})));
    Modelica.Mechanics.MultiBody.Parts.FixedRotation fixedRotation(
      animation=false,
      n(displayUnit="1") = {0,1,0},
      angle=90)
      annotation (Placement(transformation(extent={{94,-72},{114,-52}})));
    IMU IMU_sensor(start=start)
      annotation (Placement(transformation(extent={{126,-72},{164,-46}})));
    Modelica.Mechanics.MultiBody.Visualizers.FixedShape fixedShape(
      shapeType="modelica://AGILEX_BUNKER_PRO/chassis.dxf",
      length=1e-3,
      width=1e-3,
      height=1e-3,
      color={215,215,215},
      extra=1) annotation (Placement(transformation(extent={{-2,-28},{18,-8}})));
  equation
    connect(trackterrain_grouser_R.frame_b, contact_t_R.frame_b) annotation (
        Line(
        points={{137.36,4},{153,4},{153,3.76},{168.58,3.76}},
        color={95,95,95},
        thickness=0.5));
    connect(contact_t_L.frame_b, trackterrain_grouser_L.frame_b) annotation (
        Line(
        points={{-124.58,5.76},{-96,5.76},{-96,10},{-95.36,10}},
        color={95,95,95},
        thickness=0.5));
    connect(speed_L.flange, trackterrain_grouser_L.flange_a)
      annotation (Line(points={{-80,62},{-64,62},{-64,35.44}}, color={0,0,0}));
    connect(speed_R.flange, trackterrain_grouser_R.flange_a)
      annotation (Line(points={{76,64},{106,64},{106,29.44}}, color={0,0,0}));
    connect(Chassis.frame_b, trackterrain_grouser_R.frame_a) annotation (Line(
        points={{28,10},{28,14},{34,14},{34,3.76},{74,3.76}},
        color={95,95,95},
        thickness=0.5));
    connect(Chassis.frame_a, trackterrain_grouser_L.frame_a) annotation (Line(
        points={{-4,10},{-4,14},{-20,14},{-20,9.76},{-32,9.76}},
        color={95,95,95},
        thickness=0.5));
    connect(Chassis.frame_b, fixedTranslation.frame_a) annotation (Line(
        points={{28,10},{28,14},{34,14},{34,-40},{40,-40}},
        color={95,95,95},
        thickness=0.5));
    connect(combiTimeTable_L.y[1], speed_L.w_ref) annotation (Line(points={{
            -138.5,63},{-120,63},{-120,62},{-110.8,62}}, color={0,0,127}));
    connect(combiTimeTable_R.y[1], speed_R.w_ref) annotation (Line(points={{
            15.5,63},{15.5,62},{38,62},{38,64},{45.2,64}}, color={0,0,127}));
    connect(fixedTranslation.frame_b, fixedRotation.frame_a) annotation (Line(
        points={{60,-40},{88,-40},{88,-62},{94,-62}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedRotation.frame_b, IMU_sensor.frame_a) annotation (Line(
        points={{114,-62},{120,-62},{120,-58.74},{126.38,-58.74}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation.frame_b, Vicon_angles.frame_a) annotation (Line(
        points={{60,-40},{66,-40},{66,-60},{0,-60},{0,-79.5},{26,-79.5}},
        color={95,95,95},
        thickness=0.5));
    connect(Chassis.frame_a, fixedShape.frame_a) annotation (Line(
        points={{-4,10},{-12,10},{-12,-18},{-2,-18}},
        color={95,95,95},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=39,
        Interval=0.02,
        Tolerance=0.001,
        __Dymola_Algorithm="Cvode"));
  end Fig_15;

  model Fig_16
    parameter Real speed0_ref[1,1971] = -Modelica.Utilities.Streams.readRealMatrix("Fig_16_speeds.mat","speed0", 1,1971);
    parameter Real speed1_ref[1,1971] = -Modelica.Utilities.Streams.readRealMatrix("Fig_16_speeds.mat","speed1", 1,1971);
    parameter Real t[1,1971] = Modelica.Utilities.Streams.readRealMatrix("Fig_16_speeds.mat","t", 1,1971);
    parameter Real table_speed0[1971,2] = transpose([t; speed0_ref]);
    parameter Real table_speed1[1971,2] = transpose([t; speed1_ref]);
    parameter Modelica.Units.SI.Angle a = 0;
    parameter Modelica.Units.SI.Angle b = -0.0237739431485454+4.47493047320627e-06;
    parameter Modelica.Units.SI.Angle g = 0;
    final parameter Modelica.Mechanics.MultiBody.Frames.Orientation R = Modelica.Mechanics.MultiBody.Frames.axesRotations({1,2,3}, {a,b,g}, zeros(3));
    final parameter Modelica.Units.SI.Position r[3] = Modelica.Mechanics.MultiBody.Frames.resolve1(R, {0,0,0.574});
    parameter Modelica.Units.SI.Position start[3] = {3.94715324412153+0.0127822495536560, 0.45838,3.70264060014460+0.000298805176950356};
    final parameter Modelica.Units.SI.Position start1[3] = start + Modelica.Mechanics.MultiBody.Frames.resolve1(R, - {0.56507 - 0.42575 - 0.10525,0.29023,-0.574*0.5}-{0,0,0.574});
    final parameter Modelica.Units.SI.Position start2[3] = start + Modelica.Mechanics.MultiBody.Frames.resolve1(R, {-0.18554,-0.07342,0}-{0.21488,-0.04227,0}-{0.56507 - 0.42575 - 0.10525,0.29023,-0.574*0.5});
    extends Modelica.Icons.Example;
    trackterrain_grouser trackterrain_grouser_R(
      start=start2,
      a=a,
      b=b,
      g=g,
      dh={contact_t_R.contact_terrain[i].dh for i in 1:80})
      annotation (Placement(transformation(extent={{74,-8},{138,16}})));
    Modelica.Mechanics.MultiBody.Sensors.AbsoluteAngles Vicon_angles
      annotation (Placement(visible=true, transformation(
          origin={40,-79.5},
          extent={{-14,-13.5},{14,13.5}},
          rotation=0)));
    contact_t_grouser contact_t_R(
      nz=2,
      there_is_grouser={true,true,true,true,true,false,true,false,true,false,
          true,false,true,false,true,false,true,true,true,true,true,true,true,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,true,true,true,true,true,true,true,true,true,false,true,false,
          true,false,true,false,true,false,true,false,true,true,true,true,true,
          true,true,true,false,true,false,true,false,true,false,true,false,true,
          false,true,true,true,true,true},
      cm1={0.0175675,0.0117325,0.0190325,0.0131975,0.0204975,0,0.019573828,0,
          0.01634164,0,0.013109453,0,0.009877265,0,0.0066450783,0,0.00807,
          0.01537,0.009535,0.016835,0.011,0.0183,0.012465,0.019765,0,
          0.021189922,0,0.017957734,0,0.014725547,0,0.01149336,0,0.008261172,0,
          0.0073375,0.0146375,0.0088025,0.0161025,0.0102675,0.0175675,0.0117325,
          0.0190325,0.0131975,0.0204975,0,0.019573828,0,0.01634164,0,
          0.013109453,0,0.009877265,0,0.0066450783,0,0.00807,0.01537,0.009535,
          0.016835,0.011,0.0183,0.012465,0.019765,0,0.021189922,0,0.017957734,0,
          0.014725547,0,0.01149336,0,0.008261172,0,0.0073375,0.0146375,
          0.0088025,0.0161025,0.0102675},
      cm2={0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,0.01172,0,
          0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,0.0029175,
          0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,0.010255,0,
          0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,0.002185,
          0.024185,0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,
          0.01172,0,0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,
          0.0029175,0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,
          0.010255,0,0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,
          0.002185,0.024185},
      lg={0.0073,0.00437,0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,
          0.0146,0,0.0146,0,0.011695,0.002905,0.008765,0.005835,0.005835,
          0.008765,0.002905,0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,
          0.0146,0,0.01316,0.00144,0.01023,0.00437,0.0073,0.0073,0.00437,
          0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,
          0,0.011695,0.002905,0.008765,0.005835,0.005835,0.008765,0.002905,
          0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.01316,
          0.00144,0.01023,0.00437,0.0073},
      fr={true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true,false,true,false,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,true,false,true,false,false,true,false,true,
          false,true,false,true,false,true,false,true,false,true,false,true,
          false,true,false,false,true,false,true,false},
      fl={false,true,false,true,false,false,true,false,true,false,true,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          false,true,false,true,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true,false,false,true,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true})
      annotation (Placement(transformation(extent={{172,-16},{210,22}})));

    contact_t_grouser contact_t_L(
      nz=2,
      there_is_grouser={true,true,true,true,true,false,true,false,true,false,
          true,false,true,false,true,false,true,true,true,true,true,true,true,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,true,true,true,true,true,true,true,true,true,false,true,false,
          true,false,true,false,true,false,true,false,true,true,true,true,true,
          true,true,true,false,true,false,true,false,true,false,true,false,true,
          false,true,true,true,true,true},
      cm1={0.0175675,0.0117325,0.0190325,0.0131975,0.0204975,0,0.019573828,0,
          0.01634164,0,0.013109453,0,0.009877265,0,0.0066450783,0,0.00807,
          0.01537,0.009535,0.016835,0.011,0.0183,0.012465,0.019765,0,
          0.021189922,0,0.017957734,0,0.014725547,0,0.01149336,0,0.008261172,0,
          0.0073375,0.0146375,0.0088025,0.0161025,0.0102675,0.0175675,0.0117325,
          0.0190325,0.0131975,0.0204975,0,0.019573828,0,0.01634164,0,
          0.013109453,0,0.009877265,0,0.0066450783,0,0.00807,0.01537,0.009535,
          0.016835,0.011,0.0183,0.012465,0.019765,0,0.021189922,0,0.017957734,0,
          0.014725547,0,0.01149336,0,0.008261172,0,0.0073375,0.0146375,
          0.0088025,0.0161025,0.0102675},
      cm2={0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,0.01172,0,
          0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,0.0029175,
          0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,0.010255,0,
          0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,0.002185,
          0.024185,0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,
          0.01172,0,0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,
          0.0029175,0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,
          0.010255,0,0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,
          0.002185,0.024185},
      lg={0.0073,0.00437,0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,
          0.0146,0,0.0146,0,0.011695,0.002905,0.008765,0.005835,0.005835,
          0.008765,0.002905,0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,
          0.0146,0,0.01316,0.00144,0.01023,0.00437,0.0073,0.0073,0.00437,
          0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,
          0,0.011695,0.002905,0.008765,0.005835,0.005835,0.008765,0.002905,
          0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.01316,
          0.00144,0.01023,0.00437,0.0073},
      fr={true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true,false,true,false,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,true,false,true,false,false,true,false,true,
          false,true,false,true,false,true,false,true,false,true,false,true,
          false,true,false,false,true,false,true,false},
      fl={false,true,false,true,false,false,true,false,true,false,true,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          false,true,false,true,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true,false,false,true,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true})
      annotation (Placement(transformation(extent={{-128,-14},{-166,24}})));

    inner Modelica.Mechanics.MultiBody.World world
      annotation (Placement(transformation(extent={{-100,-60},{-68,-28}})));
    trackterrain_grouser trackterrain_grouser_L(
      start=start2 - r,
      a=a,
      b=b,
      g=g,
      dh={contact_t_L.contact_terrain[i].dh for i in 1:80})
      annotation (Placement(transformation(extent={{-32,-2},{-96,22}})));
    Modelica.Mechanics.Rotational.Sources.Speed speed_L
      annotation (Placement(transformation(extent={{-108,48},{-80,76}})));
    Modelica.Mechanics.Rotational.Sources.Speed speed_R
      annotation (Placement(transformation(extent={{48,50},{76,78}})));
    Modelica.Mechanics.MultiBody.Parts.BodyShape Chassis(
      animation=false,
      animateSphere=false,
      r={0,0,0.574},
      r_CM={0.03535770079,0.12598,0.574*0.5},
      m=149.4571417177305,
      I_11=(1/12)*149.4571417177305*(0.424^2 + 0.3285^2),
      I_22=(1/12)*149.4571417177305*(0.424^2 + (0.5*0.820)^2) +
          110.347597567632*0.107287700790000^2 + 39.1095441500988*
          0.302712299210000^2,
      I_33=(1/12)*149.4571417177305*((0.5*0.820)^2 + 0.3285^2) +
          110.347597567632*0.107287700790000^2 + 39.1095441500988*
          0.302712299210000^2,
      r_0(fixed=true, start=start1),
      angles_fixed=true,
      angles_start={a,b,g},
      shapeType="box",
      r_shape={0.13307,0.12598,0.075},
      lengthDirection(displayUnit="1") = {0,0,1},
      length=0.424,
      width=0.3285,
      height=0.820,
      useQuaternions=false)
      annotation (Placement(transformation(extent={{-4,-6},{28,26}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation(animation
        =false,                                                          r={
          0.56507 - 0.42575 - 0.10525,0.29023,-0.574*0.5})
      annotation (Placement(transformation(extent={{40,-50},{60,-30}})));
    inner Surface_base Surface(
      nu=2,
      nv=2,
      multiColoredSurface=false,
      color={0,255,0},
      nu1=2,
      nv1=2,
      LengthX=14.8,
      LengthZ=7.3,
      OffsetX=-2.5,
      OffsetZ=-1.7)
                 annotation (Placement(visible=true, transformation(
          origin={147,72.5},
          extent={{-15,-15.5},{15,15.5}},
          rotation=0)));
    Modelica.Blocks.Sources.CombiTimeTable combiTimeTable_L(table=table_speed1,
        timeEvents=Modelica.Blocks.Types.TimeEvents.NoTimeEvents)
      annotation (Placement(transformation(extent={{-170,48},{-140,78}})));
    Modelica.Blocks.Sources.CombiTimeTable combiTimeTable_R(table=table_speed0,
        timeEvents=Modelica.Blocks.Types.TimeEvents.NoTimeEvents)
      annotation (Placement(transformation(extent={{-16,48},{14,78}})));
    Modelica.Mechanics.MultiBody.Parts.FixedRotation fixedRotation(
      animation=false,
      n(displayUnit="1") = {0,1,0},
      angle=90)
      annotation (Placement(transformation(extent={{94,-72},{114,-52}})));
    IMU IMU_sensor(start=start)
      annotation (Placement(transformation(extent={{126,-72},{164,-46}})));
    Modelica.Mechanics.MultiBody.Visualizers.FixedShape fixedShape(
      shapeType="modelica://AGILEX_BUNKER_PRO/chassis.dxf",
      length=1e-3,
      width=1e-3,
      height=1e-3,
      color={215,215,215},
      extra=1) annotation (Placement(transformation(extent={{-2,-28},{18,-8}})));
  equation
    connect(trackterrain_grouser_R.frame_b, contact_t_R.frame_b) annotation (
        Line(
        points={{137.36,4},{153,4},{153,3.76},{168.58,3.76}},
        color={95,95,95},
        thickness=0.5));
    connect(contact_t_L.frame_b, trackterrain_grouser_L.frame_b) annotation (
        Line(
        points={{-124.58,5.76},{-96,5.76},{-96,10},{-95.36,10}},
        color={95,95,95},
        thickness=0.5));
    connect(speed_L.flange, trackterrain_grouser_L.flange_a)
      annotation (Line(points={{-80,62},{-64,62},{-64,35.44}}, color={0,0,0}));
    connect(speed_R.flange, trackterrain_grouser_R.flange_a)
      annotation (Line(points={{76,64},{106,64},{106,29.44}}, color={0,0,0}));
    connect(Chassis.frame_b, trackterrain_grouser_R.frame_a) annotation (Line(
        points={{28,10},{28,14},{34,14},{34,3.76},{74,3.76}},
        color={95,95,95},
        thickness=0.5));
    connect(Chassis.frame_a, trackterrain_grouser_L.frame_a) annotation (Line(
        points={{-4,10},{-4,14},{-20,14},{-20,9.76},{-32,9.76}},
        color={95,95,95},
        thickness=0.5));
    connect(Chassis.frame_b, fixedTranslation.frame_a) annotation (Line(
        points={{28,10},{28,14},{34,14},{34,-40},{40,-40}},
        color={95,95,95},
        thickness=0.5));
    connect(combiTimeTable_L.y[1], speed_L.w_ref) annotation (Line(points={{
            -138.5,63},{-120,63},{-120,62},{-110.8,62}}, color={0,0,127}));
    connect(combiTimeTable_R.y[1], speed_R.w_ref) annotation (Line(points={{
            15.5,63},{15.5,62},{38,62},{38,64},{45.2,64}}, color={0,0,127}));
    connect(fixedTranslation.frame_b, fixedRotation.frame_a) annotation (Line(
        points={{60,-40},{88,-40},{88,-62},{94,-62}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedRotation.frame_b, IMU_sensor.frame_a) annotation (Line(
        points={{114,-62},{120,-62},{120,-58.74},{126.38,-58.74}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation.frame_b, Vicon_angles.frame_a) annotation (Line(
        points={{60,-40},{66,-40},{66,-60},{0,-60},{0,-79.5},{26,-79.5}},
        color={95,95,95},
        thickness=0.5));
    connect(Chassis.frame_a, fixedShape.frame_a) annotation (Line(
        points={{-4,10},{-12,10},{-12,-18},{-2,-18}},
        color={95,95,95},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=39,
        Interval=0.02,
        Tolerance=0.001,
        __Dymola_Algorithm="Cvode"));
  end Fig_16;

  model Fig_18
    parameter Real speed0_ref[1,1951] = -Modelica.Utilities.Streams.readRealMatrix("Fig_18_speeds.mat","speed0", 1,1951);
    parameter Real speed1_ref[1,1951] = -Modelica.Utilities.Streams.readRealMatrix("Fig_18_speeds.mat","speed1", 1,1951);
    parameter Real t[1,1951] = Modelica.Utilities.Streams.readRealMatrix("Fig_18_speeds.mat","t", 1,1951);
    parameter Real table_speed0[1951,2] = transpose([t; speed0_ref]);
    parameter Real table_speed1[1951,2] = transpose([t; speed1_ref]);
    parameter Modelica.Units.SI.Angle a = 0;
    parameter Modelica.Units.SI.Angle b = 2.58039403953889;
    parameter Modelica.Units.SI.Angle g = 0;
    final parameter Modelica.Mechanics.MultiBody.Frames.Orientation R = Modelica.Mechanics.MultiBody.Frames.axesRotations({1,2,3}, {a,b,g}, zeros(3));
    final parameter Modelica.Units.SI.Position r[3] = Modelica.Mechanics.MultiBody.Frames.resolve1(R, {0,0,0.574});
    parameter Modelica.Units.SI.Position start[3] = {-0.00838931800339814, 0.45838, -0.00527926685977642};
    final parameter Modelica.Units.SI.Position start1[3] = start + Modelica.Mechanics.MultiBody.Frames.resolve1(R, - {0.56507 - 0.42575 - 0.10525,0.29023,-0.574*0.5}-{0,0,0.574});
    final parameter Modelica.Units.SI.Position start2[3] = start + Modelica.Mechanics.MultiBody.Frames.resolve1(R, {-0.18554,-0.07342,0}-{0.21488,-0.04227,0}-{0.56507 - 0.42575 - 0.10525,0.29023,-0.574*0.5});
    extends Modelica.Icons.Example;
    trackterrain_grouser trackterrain_grouser_R(
      start=start2,
      a=a,
      b=b,
      g=g,
      dh={contact_t_R.contact_terrain[i].dh for i in 1:80})
      annotation (Placement(transformation(extent={{74,-8},{138,16}})));
    Modelica.Mechanics.MultiBody.Sensors.AbsoluteAngles Vicon_angles
      annotation (Placement(visible=true, transformation(
          origin={40,-79.5},
          extent={{-14,-13.5},{14,13.5}},
          rotation=0)));
    contact_t_grouser contact_t_R(
      nz=2,
      there_is_grouser={true,true,true,true,true,false,true,false,true,false,
          true,false,true,false,true,false,true,true,true,true,true,true,true,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,true,true,true,true,true,true,true,true,true,false,true,false,
          true,false,true,false,true,false,true,false,true,true,true,true,true,
          true,true,true,false,true,false,true,false,true,false,true,false,true,
          false,true,true,true,true,true},
      cm1={0.0175675,0.0117325,0.0190325,0.0131975,0.0204975,0,0.019573828,0,
          0.01634164,0,0.013109453,0,0.009877265,0,0.0066450783,0,0.00807,
          0.01537,0.009535,0.016835,0.011,0.0183,0.012465,0.019765,0,
          0.021189922,0,0.017957734,0,0.014725547,0,0.01149336,0,0.008261172,0,
          0.0073375,0.0146375,0.0088025,0.0161025,0.0102675,0.0175675,0.0117325,
          0.0190325,0.0131975,0.0204975,0,0.019573828,0,0.01634164,0,
          0.013109453,0,0.009877265,0,0.0066450783,0,0.00807,0.01537,0.009535,
          0.016835,0.011,0.0183,0.012465,0.019765,0,0.021189922,0,0.017957734,0,
          0.014725547,0,0.01149336,0,0.008261172,0,0.0073375,0.0146375,
          0.0088025,0.0161025,0.0102675},
      cm2={0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,0.01172,0,
          0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,0.0029175,
          0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,0.010255,0,
          0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,0.002185,
          0.024185,0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,
          0.01172,0,0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,
          0.0029175,0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,
          0.010255,0,0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,
          0.002185,0.024185},
      lg={0.0073,0.00437,0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,
          0.0146,0,0.0146,0,0.011695,0.002905,0.008765,0.005835,0.005835,
          0.008765,0.002905,0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,
          0.0146,0,0.01316,0.00144,0.01023,0.00437,0.0073,0.0073,0.00437,
          0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,
          0,0.011695,0.002905,0.008765,0.005835,0.005835,0.008765,0.002905,
          0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.01316,
          0.00144,0.01023,0.00437,0.0073},
      fr={true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true,false,true,false,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,true,false,true,false,false,true,false,true,
          false,true,false,true,false,true,false,true,false,true,false,true,
          false,true,false,false,true,false,true,false},
      fl={false,true,false,true,false,false,true,false,true,false,true,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          false,true,false,true,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true,false,false,true,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true})
      annotation (Placement(transformation(extent={{172,-16},{210,22}})));

    contact_t_grouser contact_t_L(
      nz=2,
      there_is_grouser={true,true,true,true,true,false,true,false,true,false,
          true,false,true,false,true,false,true,true,true,true,true,true,true,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,true,true,true,true,true,true,true,true,true,false,true,false,
          true,false,true,false,true,false,true,false,true,true,true,true,true,
          true,true,true,false,true,false,true,false,true,false,true,false,true,
          false,true,true,true,true,true},
      cm1={0.0175675,0.0117325,0.0190325,0.0131975,0.0204975,0,0.019573828,0,
          0.01634164,0,0.013109453,0,0.009877265,0,0.0066450783,0,0.00807,
          0.01537,0.009535,0.016835,0.011,0.0183,0.012465,0.019765,0,
          0.021189922,0,0.017957734,0,0.014725547,0,0.01149336,0,0.008261172,0,
          0.0073375,0.0146375,0.0088025,0.0161025,0.0102675,0.0175675,0.0117325,
          0.0190325,0.0131975,0.0204975,0,0.019573828,0,0.01634164,0,
          0.013109453,0,0.009877265,0,0.0066450783,0,0.00807,0.01537,0.009535,
          0.016835,0.011,0.0183,0.012465,0.019765,0,0.021189922,0,0.017957734,0,
          0.014725547,0,0.01149336,0,0.008261172,0,0.0073375,0.0146375,
          0.0088025,0.0161025,0.0102675},
      cm2={0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,0.01172,0,
          0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,0.0029175,
          0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,0.010255,0,
          0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,0.002185,
          0.024185,0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,
          0.01172,0,0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,
          0.0029175,0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,
          0.010255,0,0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,
          0.002185,0.024185},
      lg={0.0073,0.00437,0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,
          0.0146,0,0.0146,0,0.011695,0.002905,0.008765,0.005835,0.005835,
          0.008765,0.002905,0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,
          0.0146,0,0.01316,0.00144,0.01023,0.00437,0.0073,0.0073,0.00437,
          0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,
          0,0.011695,0.002905,0.008765,0.005835,0.005835,0.008765,0.002905,
          0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.01316,
          0.00144,0.01023,0.00437,0.0073},
      fr={true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true,false,true,false,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,true,false,true,false,false,true,false,true,
          false,true,false,true,false,true,false,true,false,true,false,true,
          false,true,false,false,true,false,true,false},
      fl={false,true,false,true,false,false,true,false,true,false,true,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          false,true,false,true,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true,false,false,true,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true})
      annotation (Placement(transformation(extent={{-128,-14},{-166,24}})));

    inner Modelica.Mechanics.MultiBody.World world
      annotation (Placement(transformation(extent={{-100,-60},{-68,-28}})));
    trackterrain_grouser trackterrain_grouser_L(
      start=start2 - r,
      a=a,
      b=b,
      g=g,
      dh={contact_t_L.contact_terrain[i].dh for i in 1:80})
      annotation (Placement(transformation(extent={{-32,-2},{-96,22}})));
    Modelica.Mechanics.Rotational.Sources.Speed speed_L
      annotation (Placement(transformation(extent={{-108,48},{-80,76}})));
    Modelica.Mechanics.Rotational.Sources.Speed speed_R
      annotation (Placement(transformation(extent={{48,50},{76,78}})));
    Modelica.Mechanics.MultiBody.Parts.BodyShape Chassis(
      animation=false,
      animateSphere=false,
      r={0,0,0.574},
      r_CM={0.03535770079,0.12598,0.574*0.5},
      m=149.4571417177305,
      I_11=(1/12)*149.4571417177305*(0.424^2 + 0.3285^2),
      I_22=(1/12)*149.4571417177305*(0.424^2 + (0.5*0.820)^2) +
          110.347597567632*0.107287700790000^2 + 39.1095441500988*
          0.302712299210000^2,
      I_33=(1/12)*149.4571417177305*((0.5*0.820)^2 + 0.3285^2) +
          110.347597567632*0.107287700790000^2 + 39.1095441500988*
          0.302712299210000^2,
      r_0(fixed=true, start=start1),
      angles_fixed=true,
      angles_start={a,b,g},
      shapeType="box",
      r_shape={0.13307,0.12598,0.075},
      lengthDirection(displayUnit="1") = {0,0,1},
      length=0.424,
      width=0.3285,
      height=0.820,
      useQuaternions=false)
      annotation (Placement(transformation(extent={{-4,-6},{28,26}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation(animation
        =false,                                                          r={
          0.56507 - 0.42575 - 0.10525,0.29023,-0.574*0.5})
      annotation (Placement(transformation(extent={{40,-50},{60,-30}})));
    inner Surface_base Surface(
      nu=2,
      nv=2,
      multiColoredSurface=false,
      color={0,255,0},
      DefaultCharSurface=false,
      nu1=2,
      nv1=2,
      LengthX=14.8,
      LengthZ=7.2,
      mu=[0.9],
      kc=[52530],
      kf=[1127970],
      kv=[530100.0],
      n=[0.9],
      c=[0],
      tanf=[tan((20/180)*Modelica.Constants.pi)],
      K=[0.025],
      gammas=[11000])
                 annotation (Placement(visible=true, transformation(
          origin={147,72.5},
          extent={{-15,-15.5},{15,15.5}},
          rotation=0)));
    Modelica.Blocks.Sources.CombiTimeTable combiTimeTable_L(table=table_speed1,
        timeEvents=Modelica.Blocks.Types.TimeEvents.NoTimeEvents)
      annotation (Placement(transformation(extent={{-170,48},{-140,78}})));
    Modelica.Blocks.Sources.CombiTimeTable combiTimeTable_R(table=table_speed0,
        timeEvents=Modelica.Blocks.Types.TimeEvents.NoTimeEvents)
      annotation (Placement(transformation(extent={{-16,48},{14,78}})));
    Modelica.Mechanics.MultiBody.Parts.FixedRotation fixedRotation(
      animation=false,
      n(displayUnit="1") = {0,1,0},
      angle=90)
      annotation (Placement(transformation(extent={{94,-72},{114,-52}})));
    IMU IMU_sensor(start=start)
      annotation (Placement(transformation(extent={{126,-72},{164,-46}})));
    Modelica.Mechanics.MultiBody.Visualizers.FixedShape fixedShape(
      shapeType="modelica://AGILEX_BUNKER_PRO/chassis.dxf",
      length=1e-3,
      width=1e-3,
      height=1e-3,
      color={215,215,215},
      extra=1) annotation (Placement(transformation(extent={{-2,-28},{18,-8}})));
  equation
    connect(trackterrain_grouser_R.frame_b, contact_t_R.frame_b) annotation (
        Line(
        points={{137.36,4},{153,4},{153,3.76},{168.58,3.76}},
        color={95,95,95},
        thickness=0.5));
    connect(contact_t_L.frame_b, trackterrain_grouser_L.frame_b) annotation (
        Line(
        points={{-124.58,5.76},{-96,5.76},{-96,10},{-95.36,10}},
        color={95,95,95},
        thickness=0.5));
    connect(speed_L.flange, trackterrain_grouser_L.flange_a)
      annotation (Line(points={{-80,62},{-64,62},{-64,35.44}}, color={0,0,0}));
    connect(speed_R.flange, trackterrain_grouser_R.flange_a)
      annotation (Line(points={{76,64},{106,64},{106,29.44}}, color={0,0,0}));
    connect(Chassis.frame_b, trackterrain_grouser_R.frame_a) annotation (Line(
        points={{28,10},{28,14},{34,14},{34,3.76},{74,3.76}},
        color={95,95,95},
        thickness=0.5));
    connect(Chassis.frame_a, trackterrain_grouser_L.frame_a) annotation (Line(
        points={{-4,10},{-4,14},{-20,14},{-20,9.76},{-32,9.76}},
        color={95,95,95},
        thickness=0.5));
    connect(Chassis.frame_b, fixedTranslation.frame_a) annotation (Line(
        points={{28,10},{28,14},{34,14},{34,-40},{40,-40}},
        color={95,95,95},
        thickness=0.5));
    connect(combiTimeTable_L.y[1], speed_L.w_ref) annotation (Line(points={{
            -138.5,63},{-120,63},{-120,62},{-110.8,62}}, color={0,0,127}));
    connect(combiTimeTable_R.y[1], speed_R.w_ref) annotation (Line(points={{
            15.5,63},{15.5,62},{38,62},{38,64},{45.2,64}}, color={0,0,127}));
    connect(fixedTranslation.frame_b, fixedRotation.frame_a) annotation (Line(
        points={{60,-40},{88,-40},{88,-62},{94,-62}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedRotation.frame_b, IMU_sensor.frame_a) annotation (Line(
        points={{114,-62},{120,-62},{120,-58.74},{126.38,-58.74}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation.frame_b, Vicon_angles.frame_a) annotation (Line(
        points={{60,-40},{66,-40},{66,-60},{0,-60},{0,-79.5},{26,-79.5}},
        color={95,95,95},
        thickness=0.5));
    connect(Chassis.frame_a, fixedShape.frame_a) annotation (Line(
        points={{-4,10},{-12,10},{-12,-18},{-2,-18}},
        color={95,95,95},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=38.5,
        Interval=0.02,
        Tolerance=0.001,
        __Dymola_Algorithm="Cvode"));
  end Fig_18;

  model Fig_20
    parameter Real speed0_ref[1,1901] = -Modelica.Utilities.Streams.readRealMatrix("Fig_20_speeds.mat","speed0", 1,1901);
    parameter Real speed1_ref[1,1901] = -Modelica.Utilities.Streams.readRealMatrix("Fig_20_speeds.mat","speed1", 1,1901);
    parameter Real t[1,1901] = Modelica.Utilities.Streams.readRealMatrix("Fig_20_speeds.mat","t", 1,1901);
    parameter Real table_speed0[1901,2] = transpose([t; speed0_ref]);
    parameter Real table_speed1[1901,2] = transpose([t; speed1_ref]);
    parameter Modelica.Units.SI.Angle a = 0;
    parameter Modelica.Units.SI.Angle b = 2.73916727742618-0.1;
    parameter Modelica.Units.SI.Angle g = 0;
    final parameter Modelica.Mechanics.MultiBody.Frames.Orientation R = Modelica.Mechanics.MultiBody.Frames.axesRotations({1,2,3}, {a,b,g}, zeros(3));
    final parameter Modelica.Units.SI.Position r[3] = Modelica.Mechanics.MultiBody.Frames.resolve1(R, {0,0,0.574});
    parameter Modelica.Units.SI.Position start[3] = {-0.00903874225186615 + 0.000928943489953018, 0.45838, -0.00496296610164084 + 0.000491427617430023};
    final parameter Modelica.Units.SI.Position start1[3] = start + Modelica.Mechanics.MultiBody.Frames.resolve1(R, - {0.56507 - 0.42575 - 0.10525,0.29023,-0.574*0.5}-{0,0,0.574});
    final parameter Modelica.Units.SI.Position start2[3] = start + Modelica.Mechanics.MultiBody.Frames.resolve1(R, {-0.18554,-0.07342,0}-{0.21488,-0.04227,0}-{0.56507 - 0.42575 - 0.10525,0.29023,-0.574*0.5});
    extends Modelica.Icons.Example;
    trackterrain_grouser trackterrain_grouser_R(
      start=start2,
      a=a,
      b=b,
      g=g,
      dh={contact_t_R.contact_terrain[i].dh for i in 1:80})
      annotation (Placement(transformation(extent={{74,-8},{138,16}})));
    Modelica.Mechanics.MultiBody.Sensors.AbsoluteAngles Vicon_angles
      annotation (Placement(visible=true, transformation(
          origin={40,-79.5},
          extent={{-14,-13.5},{14,13.5}},
          rotation=0)));
    contact_t_grouser contact_t_R(
      nz=2,
      there_is_grouser={true,true,true,true,true,false,true,false,true,false,
          true,false,true,false,true,false,true,true,true,true,true,true,true,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,true,true,true,true,true,true,true,true,true,false,true,false,
          true,false,true,false,true,false,true,false,true,true,true,true,true,
          true,true,true,false,true,false,true,false,true,false,true,false,true,
          false,true,true,true,true,true},
      cm1={0.0175675,0.0117325,0.0190325,0.0131975,0.0204975,0,0.019573828,0,
          0.01634164,0,0.013109453,0,0.009877265,0,0.0066450783,0,0.00807,
          0.01537,0.009535,0.016835,0.011,0.0183,0.012465,0.019765,0,
          0.021189922,0,0.017957734,0,0.014725547,0,0.01149336,0,0.008261172,0,
          0.0073375,0.0146375,0.0088025,0.0161025,0.0102675,0.0175675,0.0117325,
          0.0190325,0.0131975,0.0204975,0,0.019573828,0,0.01634164,0,
          0.013109453,0,0.009877265,0,0.0066450783,0,0.00807,0.01537,0.009535,
          0.016835,0.011,0.0183,0.012465,0.019765,0,0.021189922,0,0.017957734,0,
          0.014725547,0,0.01149336,0,0.008261172,0,0.0073375,0.0146375,
          0.0088025,0.0161025,0.0102675},
      cm2={0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,0.01172,0,
          0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,0.0029175,
          0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,0.010255,0,
          0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,0.002185,
          0.024185,0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,
          0.01172,0,0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,
          0.0029175,0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,
          0.010255,0,0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,
          0.002185,0.024185},
      lg={0.0073,0.00437,0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,
          0.0146,0,0.0146,0,0.011695,0.002905,0.008765,0.005835,0.005835,
          0.008765,0.002905,0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,
          0.0146,0,0.01316,0.00144,0.01023,0.00437,0.0073,0.0073,0.00437,
          0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,
          0,0.011695,0.002905,0.008765,0.005835,0.005835,0.008765,0.002905,
          0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.01316,
          0.00144,0.01023,0.00437,0.0073},
      fr={true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true,false,true,false,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,true,false,true,false,false,true,false,true,
          false,true,false,true,false,true,false,true,false,true,false,true,
          false,true,false,false,true,false,true,false},
      fl={false,true,false,true,false,false,true,false,true,false,true,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          false,true,false,true,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true,false,false,true,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true})
      annotation (Placement(transformation(extent={{172,-16},{210,22}})));

    contact_t_grouser contact_t_L(
      nz=2,
      there_is_grouser={true,true,true,true,true,false,true,false,true,false,
          true,false,true,false,true,false,true,true,true,true,true,true,true,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,true,true,true,true,true,true,true,true,true,false,true,false,
          true,false,true,false,true,false,true,false,true,true,true,true,true,
          true,true,true,false,true,false,true,false,true,false,true,false,true,
          false,true,true,true,true,true},
      cm1={0.0175675,0.0117325,0.0190325,0.0131975,0.0204975,0,0.019573828,0,
          0.01634164,0,0.013109453,0,0.009877265,0,0.0066450783,0,0.00807,
          0.01537,0.009535,0.016835,0.011,0.0183,0.012465,0.019765,0,
          0.021189922,0,0.017957734,0,0.014725547,0,0.01149336,0,0.008261172,0,
          0.0073375,0.0146375,0.0088025,0.0161025,0.0102675,0.0175675,0.0117325,
          0.0190325,0.0131975,0.0204975,0,0.019573828,0,0.01634164,0,
          0.013109453,0,0.009877265,0,0.0066450783,0,0.00807,0.01537,0.009535,
          0.016835,0.011,0.0183,0.012465,0.019765,0,0.021189922,0,0.017957734,0,
          0.014725547,0,0.01149336,0,0.008261172,0,0.0073375,0.0146375,
          0.0088025,0.0161025,0.0102675},
      cm2={0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,0.01172,0,
          0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,0.0029175,
          0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,0.010255,0,
          0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,0.002185,
          0.024185,0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,
          0.01172,0,0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,
          0.0029175,0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,
          0.010255,0,0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,
          0.002185,0.024185},
      lg={0.0073,0.00437,0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,
          0.0146,0,0.0146,0,0.011695,0.002905,0.008765,0.005835,0.005835,
          0.008765,0.002905,0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,
          0.0146,0,0.01316,0.00144,0.01023,0.00437,0.0073,0.0073,0.00437,
          0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,
          0,0.011695,0.002905,0.008765,0.005835,0.005835,0.008765,0.002905,
          0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.01316,
          0.00144,0.01023,0.00437,0.0073},
      fr={true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true,false,true,false,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,true,false,true,false,false,true,false,true,
          false,true,false,true,false,true,false,true,false,true,false,true,
          false,true,false,false,true,false,true,false},
      fl={false,true,false,true,false,false,true,false,true,false,true,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          false,true,false,true,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true,false,false,true,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true})
      annotation (Placement(transformation(extent={{-128,-14},{-166,24}})));

    inner Modelica.Mechanics.MultiBody.World world
      annotation (Placement(transformation(extent={{-100,-60},{-68,-28}})));
    trackterrain_grouser trackterrain_grouser_L(
      start=start2 - r,
      a=a,
      b=b,
      g=g,
      dh={contact_t_L.contact_terrain[i].dh for i in 1:80})
      annotation (Placement(transformation(extent={{-32,-2},{-96,22}})));
    Modelica.Mechanics.Rotational.Sources.Speed speed_L
      annotation (Placement(transformation(extent={{-108,48},{-80,76}})));
    Modelica.Mechanics.Rotational.Sources.Speed speed_R
      annotation (Placement(transformation(extent={{48,50},{76,78}})));
    Modelica.Mechanics.MultiBody.Parts.BodyShape Chassis(
      animation=false,
      animateSphere=false,
      r={0,0,0.574},
      r_CM={0.03535770079,0.12598,0.574*0.5},
      m=149.4571417177305,
      I_11=(1/12)*149.4571417177305*(0.424^2 + 0.3285^2),
      I_22=(1/12)*149.4571417177305*(0.424^2 + (0.5*0.820)^2) +
          110.347597567632*0.107287700790000^2 + 39.1095441500988*
          0.302712299210000^2,
      I_33=(1/12)*149.4571417177305*((0.5*0.820)^2 + 0.3285^2) +
          110.347597567632*0.107287700790000^2 + 39.1095441500988*
          0.302712299210000^2,
      r_0(fixed=true, start=start1),
      angles_fixed=true,
      angles_start={a,b,g},
      shapeType="box",
      r_shape={0.13307,0.12598,0.075},
      lengthDirection(displayUnit="1") = {0,0,1},
      length=0.424,
      width=0.3285,
      height=0.820,
      useQuaternions=false)
      annotation (Placement(transformation(extent={{-4,-6},{28,26}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation(animation
        =false,                                                          r={
          0.56507 - 0.42575 - 0.10525,0.29023,-0.574*0.5})
      annotation (Placement(transformation(extent={{40,-50},{60,-30}})));
    inner Surface_base Surface(
      nu=2,
      nv=2,
      multiColoredSurface=false,
      color={0,255,0},
      DefaultCharSurface=false,
      nu1=2,
      nv1=2,
      LengthX=14.8,
      LengthZ=7.2,
      mu=[0.9],
      kc=[52530],
      kf=[1127970],
      kv=[530100.0],
      n=[0.9],
      c=[0],
      tanf=[tan((20/180)*Modelica.Constants.pi)],
      K=[0.025],
      gammas=[11000])
                 annotation (Placement(visible=true, transformation(
          origin={147,72.5},
          extent={{-15,-15.5},{15,15.5}},
          rotation=0)));
    Modelica.Blocks.Sources.CombiTimeTable combiTimeTable_L(table=table_speed1,
        timeEvents=Modelica.Blocks.Types.TimeEvents.NoTimeEvents)
      annotation (Placement(transformation(extent={{-170,48},{-140,78}})));
    Modelica.Blocks.Sources.CombiTimeTable combiTimeTable_R(table=table_speed0,
        timeEvents=Modelica.Blocks.Types.TimeEvents.NoTimeEvents)
      annotation (Placement(transformation(extent={{-16,48},{14,78}})));
    Modelica.Mechanics.MultiBody.Parts.FixedRotation fixedRotation(
      animation=false,
      n(displayUnit="1") = {0,1,0},
      angle=90)
      annotation (Placement(transformation(extent={{94,-72},{114,-52}})));
    IMU IMU_sensor(start=start)
      annotation (Placement(transformation(extent={{126,-72},{164,-46}})));
    Modelica.Mechanics.MultiBody.Visualizers.FixedShape fixedShape(
      shapeType="modelica://AGILEX_BUNKER_PRO/chassis.dxf",
      length=1e-3,
      width=1e-3,
      height=1e-3,
      color={215,215,215},
      extra=1) annotation (Placement(transformation(extent={{-2,-28},{18,-8}})));
  equation
    connect(trackterrain_grouser_R.frame_b, contact_t_R.frame_b) annotation (
        Line(
        points={{137.36,4},{153,4},{153,3.76},{168.58,3.76}},
        color={95,95,95},
        thickness=0.5));
    connect(contact_t_L.frame_b, trackterrain_grouser_L.frame_b) annotation (
        Line(
        points={{-124.58,5.76},{-96,5.76},{-96,10},{-95.36,10}},
        color={95,95,95},
        thickness=0.5));
    connect(speed_L.flange, trackterrain_grouser_L.flange_a)
      annotation (Line(points={{-80,62},{-64,62},{-64,35.44}}, color={0,0,0}));
    connect(speed_R.flange, trackterrain_grouser_R.flange_a)
      annotation (Line(points={{76,64},{106,64},{106,29.44}}, color={0,0,0}));
    connect(Chassis.frame_b, trackterrain_grouser_R.frame_a) annotation (Line(
        points={{28,10},{28,14},{34,14},{34,3.76},{74,3.76}},
        color={95,95,95},
        thickness=0.5));
    connect(Chassis.frame_a, trackterrain_grouser_L.frame_a) annotation (Line(
        points={{-4,10},{-4,14},{-20,14},{-20,9.76},{-32,9.76}},
        color={95,95,95},
        thickness=0.5));
    connect(Chassis.frame_b, fixedTranslation.frame_a) annotation (Line(
        points={{28,10},{28,14},{34,14},{34,-40},{40,-40}},
        color={95,95,95},
        thickness=0.5));
    connect(combiTimeTable_L.y[1], speed_L.w_ref) annotation (Line(points={{
            -138.5,63},{-120,63},{-120,62},{-110.8,62}}, color={0,0,127}));
    connect(combiTimeTable_R.y[1], speed_R.w_ref) annotation (Line(points={{
            15.5,63},{15.5,62},{38,62},{38,64},{45.2,64}}, color={0,0,127}));
    connect(fixedTranslation.frame_b, fixedRotation.frame_a) annotation (Line(
        points={{60,-40},{88,-40},{88,-62},{94,-62}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedRotation.frame_b, IMU_sensor.frame_a) annotation (Line(
        points={{114,-62},{120,-62},{120,-58.74},{126.38,-58.74}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation.frame_b, Vicon_angles.frame_a) annotation (Line(
        points={{60,-40},{66,-40},{66,-60},{0,-60},{0,-79.5},{26,-79.5}},
        color={95,95,95},
        thickness=0.5));
    connect(Chassis.frame_a, fixedShape.frame_a) annotation (Line(
        points={{-4,10},{-12,10},{-12,-18},{-2,-18}},
        color={95,95,95},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=37.5,
        Interval=0.02,
        Tolerance=0.001,
        __Dymola_Algorithm="Cvode"));
  end Fig_20;

  model Fig_21
    parameter Modelica.Units.SI.Angle a = 0;
    parameter Modelica.Units.SI.Angle b = 0;
    parameter Modelica.Units.SI.Angle g = 0;
    final parameter Modelica.Mechanics.MultiBody.Frames.Orientation R = Modelica.Mechanics.MultiBody.Frames.axesRotations({1,2,3}, {a,b,g}, zeros(3));
    final parameter Modelica.Units.SI.Position r[3] = Modelica.Mechanics.MultiBody.Frames.resolve1(R, {0,0,0.574});
    parameter Modelica.Units.SI.Position start[3] = {-1.0601, 0.45838, 0};
    final parameter Modelica.Units.SI.Position start1[3] = start + Modelica.Mechanics.MultiBody.Frames.resolve1(R, - {0.56507 - 0.42575 - 0.10525,0.29023,-0.574*0.5}-{0,0,0.574});
    final parameter Modelica.Units.SI.Position start2[3] = start + Modelica.Mechanics.MultiBody.Frames.resolve1(R, {-0.18554,-0.07342,0}-{0.21488,-0.04227,0}-{0.56507 - 0.42575 - 0.10525,0.29023,-0.574*0.5});
    extends Modelica.Icons.Example;
    trackterrain_grouser trackterrain_grouser_R(
      start=start2,
      a=a,
      b=b,
      g=g,
      dh={contact_t_R.contact_terrain[i].dh for i in 1:80})
      annotation (Placement(transformation(extent={{74,-8},{138,16}})));
    Modelica.Mechanics.MultiBody.Sensors.AbsoluteAngles Vicon_angles
      annotation (Placement(visible=true, transformation(
          origin={40,-79.5},
          extent={{-14,-13.5},{14,13.5}},
          rotation=0)));
    contact_t_grouser contact_t_R(
      nz=2,
      there_is_grouser={true,true,true,true,true,false,true,false,true,false,
          true,false,true,false,true,false,true,true,true,true,true,true,true,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,true,true,true,true,true,true,true,true,true,false,true,false,
          true,false,true,false,true,false,true,false,true,true,true,true,true,
          true,true,true,false,true,false,true,false,true,false,true,false,true,
          false,true,true,true,true,true},
      cm1={0.0175675,0.0117325,0.0190325,0.0131975,0.0204975,0,0.019573828,0,
          0.01634164,0,0.013109453,0,0.009877265,0,0.0066450783,0,0.00807,
          0.01537,0.009535,0.016835,0.011,0.0183,0.012465,0.019765,0,
          0.021189922,0,0.017957734,0,0.014725547,0,0.01149336,0,0.008261172,0,
          0.0073375,0.0146375,0.0088025,0.0161025,0.0102675,0.0175675,0.0117325,
          0.0190325,0.0131975,0.0204975,0,0.019573828,0,0.01634164,0,
          0.013109453,0,0.009877265,0,0.0066450783,0,0.00807,0.01537,0.009535,
          0.016835,0.011,0.0183,0.012465,0.019765,0,0.021189922,0,0.017957734,0,
          0.014725547,0,0.01149336,0,0.008261172,0,0.0073375,0.0146375,
          0.0088025,0.0161025,0.0102675},
      cm2={0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,0.01172,0,
          0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,0.0029175,
          0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,0.010255,0,
          0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,0.002185,
          0.024185,0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,
          0.01172,0,0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,
          0.0029175,0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,
          0.010255,0,0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,
          0.002185,0.024185},
      lg={0.0073,0.00437,0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,
          0.0146,0,0.0146,0,0.011695,0.002905,0.008765,0.005835,0.005835,
          0.008765,0.002905,0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,
          0.0146,0,0.01316,0.00144,0.01023,0.00437,0.0073,0.0073,0.00437,
          0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,
          0,0.011695,0.002905,0.008765,0.005835,0.005835,0.008765,0.002905,
          0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.01316,
          0.00144,0.01023,0.00437,0.0073},
      fr={true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true,false,true,false,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,true,false,true,false,false,true,false,true,
          false,true,false,true,false,true,false,true,false,true,false,true,
          false,true,false,false,true,false,true,false},
      fl={false,true,false,true,false,false,true,false,true,false,true,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          false,true,false,true,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true,false,false,true,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true})
      annotation (Placement(transformation(extent={{172,-16},{210,22}})));

    contact_t_grouser contact_t_L(
      nz=2,
      there_is_grouser={true,true,true,true,true,false,true,false,true,false,
          true,false,true,false,true,false,true,true,true,true,true,true,true,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,true,true,true,true,true,true,true,true,true,false,true,false,
          true,false,true,false,true,false,true,false,true,true,true,true,true,
          true,true,true,false,true,false,true,false,true,false,true,false,true,
          false,true,true,true,true,true},
      cm1={0.0175675,0.0117325,0.0190325,0.0131975,0.0204975,0,0.019573828,0,
          0.01634164,0,0.013109453,0,0.009877265,0,0.0066450783,0,0.00807,
          0.01537,0.009535,0.016835,0.011,0.0183,0.012465,0.019765,0,
          0.021189922,0,0.017957734,0,0.014725547,0,0.01149336,0,0.008261172,0,
          0.0073375,0.0146375,0.0088025,0.0161025,0.0102675,0.0175675,0.0117325,
          0.0190325,0.0131975,0.0204975,0,0.019573828,0,0.01634164,0,
          0.013109453,0,0.009877265,0,0.0066450783,0,0.00807,0.01537,0.009535,
          0.016835,0.011,0.0183,0.012465,0.019765,0,0.021189922,0,0.017957734,0,
          0.014725547,0,0.01149336,0,0.008261172,0,0.0073375,0.0146375,
          0.0088025,0.0161025,0.0102675},
      cm2={0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,0.01172,0,
          0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,0.0029175,
          0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,0.010255,0,
          0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,0.002185,
          0.024185,0.00365,0.02565,0.005115,0.027115,0.00658,0,0.00879,0,
          0.01172,0,0.01465,0,0.01758,0,0.02051,0,0.0219875,0.0014525,0.0234525,
          0.0029175,0.0249175,0.0043825,0.0263825,0.0058475,0,0.007325,0,
          0.010255,0,0.013185,0,0.016115,0,0.019045,0,0.021255,0.00072,0.02272,
          0.002185,0.024185},
      lg={0.0073,0.00437,0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,
          0.0146,0,0.0146,0,0.011695,0.002905,0.008765,0.005835,0.005835,
          0.008765,0.002905,0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,
          0.0146,0,0.01316,0.00144,0.01023,0.00437,0.0073,0.0073,0.00437,
          0.01023,0.00144,0.01316,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,
          0,0.011695,0.002905,0.008765,0.005835,0.005835,0.008765,0.002905,
          0.011695,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.0146,0,0.01316,
          0.00144,0.01023,0.00437,0.0073},
      fr={true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true,false,true,false,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,true,false,true,false,false,true,false,true,
          false,true,false,true,false,true,false,true,false,true,false,true,
          false,true,false,false,true,false,true,false},
      fl={false,true,false,true,false,false,true,false,true,false,true,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          false,true,false,true,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true,false,false,true,false,
          true,false,true,false,true,false,true,false,true,false,true,false,
          true,false,true,false,false,true,false,true,false,true,false,true,
          false,true,false,true,false,true,false,true})
      annotation (Placement(transformation(extent={{-128,-14},{-166,24}})));

    inner Modelica.Mechanics.MultiBody.World world
      annotation (Placement(transformation(extent={{-100,-60},{-68,-28}})));
    trackterrain_grouser trackterrain_grouser_L(
      start=start2 - r,
      a=a,
      b=b,
      g=g,
      dh={contact_t_L.contact_terrain[i].dh for i in 1:80})
      annotation (Placement(transformation(extent={{-32,-2},{-96,22}})));
    Modelica.Mechanics.Rotational.Sources.Speed speed_L
      annotation (Placement(transformation(extent={{-108,48},{-80,76}})));
    Modelica.Mechanics.Rotational.Sources.Speed speed_R
      annotation (Placement(transformation(extent={{48,50},{76,78}})));
    Modelica.Mechanics.MultiBody.Parts.BodyShape Chassis(
      animation=false,
      animateSphere=false,
      r={0,0,0.574},
      r_CM={0.03535770079,0.12598,0.574*0.5},
      m=149.4571417177305,
      I_11=(1/12)*149.4571417177305*(0.424^2 + 0.3285^2),
      I_22=(1/12)*149.4571417177305*(0.424^2 + (0.5*0.820)^2) + 110.347597567632*0.107287700790000
          ^2 + 39.1095441500988*0.302712299210000^2,
      I_33=(1/12)*149.4571417177305*((0.5*0.820)^2 + 0.3285^2) + 110.347597567632*
          0.107287700790000^2 + 39.1095441500988*0.302712299210000^2,
      r_0(fixed=true, start=start1),
      angles_fixed=true,
      angles_start={a,b,g},
      shapeType="box",
      r_shape={0.13307,0.12598,0.075},
      lengthDirection(displayUnit="1") = {0,0,1},
      length=0.424,
      width=0.3285,
      height=0.820,
      useQuaternions=false)
      annotation (Placement(transformation(extent={{-4,-6},{28,26}})));
    Modelica.Mechanics.MultiBody.Parts.FixedTranslation fixedTranslation(animation
        =false,                                                          r={
          0.56507 - 0.42575 - 0.10525,0.29023,-0.574*0.5})
      annotation (Placement(transformation(extent={{40,-50},{60,-30}})));
    Modelica.Mechanics.MultiBody.Parts.FixedRotation fixedRotation(
      animation=false,
      n(displayUnit="1") = {0,1,0},
      angle=90)
      annotation (Placement(transformation(extent={{94,-72},{114,-52}})));
    IMU IMU_sensor(start=start)
      annotation (Placement(transformation(extent={{126,-72},{164,-46}})));
    inner Surface_with_step Surface(
      nu=10,
      multiColoredSurface=true,
      color={0,255,0},
      FlatSurface=false,
      DefaultCharSurface=false,
      nu1=2,
      nv1=2,
      LengthX=8,
      LengthZ=8,
      mu=[0.9])
      annotation (Placement(transformation(extent={{138,56},{158,76}})));
    Modelica.Blocks.Sources.Ramp ramp_R(
      duration=0,
      height=-5.72846083117084,
      startTime=3) annotation (Placement(visible=true, transformation(
          origin={14,63},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    Modelica.Blocks.Sources.Ramp ramp_L(
      duration=0,
      height=-5.72846083117084,
      startTime=3) annotation (Placement(visible=true, transformation(
          origin={-146,63},
          extent={{-10,-10},{10,10}},
          rotation=0)));
    Modelica.Mechanics.MultiBody.Visualizers.FixedShape fixedShape(
      shapeType="modelica://AGILEX_BUNKER_PRO/chassis.dxf",
      length=1e-3,
      width=1e-3,
      height=1e-3,
      color={215,215,215},
      extra=1) annotation (Placement(transformation(extent={{-4,-32},{16,-12}})));
  equation
    connect(trackterrain_grouser_R.frame_b, contact_t_R.frame_b) annotation (
        Line(
        points={{137.36,4},{153,4},{153,3.76},{168.58,3.76}},
        color={95,95,95},
        thickness=0.5));
    connect(contact_t_L.frame_b, trackterrain_grouser_L.frame_b) annotation (
        Line(
        points={{-124.58,5.76},{-96,5.76},{-96,10},{-95.36,10}},
        color={95,95,95},
        thickness=0.5));
    connect(speed_L.flange, trackterrain_grouser_L.flange_a)
      annotation (Line(points={{-80,62},{-64,62},{-64,35.44}}, color={0,0,0}));
    connect(speed_R.flange, trackterrain_grouser_R.flange_a)
      annotation (Line(points={{76,64},{106,64},{106,29.44}}, color={0,0,0}));
    connect(Chassis.frame_b, trackterrain_grouser_R.frame_a) annotation (Line(
        points={{28,10},{28,14},{34,14},{34,3.76},{74,3.76}},
        color={95,95,95},
        thickness=0.5));
    connect(Chassis.frame_a, trackterrain_grouser_L.frame_a) annotation (Line(
        points={{-4,10},{-4,14},{-20,14},{-20,9.76},{-32,9.76}},
        color={95,95,95},
        thickness=0.5));
    connect(Chassis.frame_b, fixedTranslation.frame_a) annotation (Line(
        points={{28,10},{28,14},{34,14},{34,-40},{40,-40}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation.frame_b, fixedRotation.frame_a) annotation (Line(
        points={{60,-40},{88,-40},{88,-62},{94,-62}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedRotation.frame_b, IMU_sensor.frame_a) annotation (Line(
        points={{114,-62},{120,-62},{120,-58.74},{126.38,-58.74}},
        color={95,95,95},
        thickness=0.5));
    connect(fixedTranslation.frame_b, Vicon_angles.frame_a) annotation (Line(
        points={{60,-40},{66,-40},{66,-60},{0,-60},{0,-79.5},{26,-79.5}},
        color={95,95,95},
        thickness=0.5));
    connect(ramp_R.y, speed_R.w_ref)
      annotation (Line(points={{25,63},{28,64},{45.2,64}}, color={0,0,127}));
    connect(ramp_L.y, speed_L.w_ref) annotation (Line(points={{-135,63},{-118,
            63},{-118,62},{-110.8,62}}, color={0,0,127}));
    connect(Chassis.frame_a, fixedShape.frame_a) annotation (Line(
        points={{-4,10},{-14,10},{-14,-22},{-4,-22}},
        color={95,95,95},
        thickness=0.5));
    annotation (Icon(coordinateSystem(preserveAspectRatio=false)), Diagram(
          coordinateSystem(preserveAspectRatio=false)),
      experiment(
        StopTime=12,
        Interval=0.02,
        Tolerance=0.001,
        __Dymola_Algorithm="Cvode"));
  end Fig_21;
  annotation (
    Icon(coordinateSystem(grid = {2, 0})),
    uses(Modelica(version = "4.0.0")));
end AGILEX_BUNKER_PRO;
