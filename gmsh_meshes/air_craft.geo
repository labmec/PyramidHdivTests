
// ---- Gmsh Script ----
// ---- Vertical wellbore wih hybrid 3D mesh  ----
// Created 01/02/2018 by Omar Durán
// Labmec, University of Campinas
// --------------------------------------------

SetFactory("OpenCASCADE");

Mesh.Algorithm=1; // 2D mesh algorithm  (1) MeshAdapt (default)(5) Delaunay (6) Frontal
Mesh.ElementOrder=1; // 1=linear elements, N (<6) = elements of higher order
Geometry.OCCSewFaces = 1;
Merge "air_craft.msh";
//Merge "benchmarks_3d_Falcon_SurfaceMeshAnisoCurvature.stl";
CreateTopology;



// points

Point(2)={-1.37177,0.592146,0.272342};
Point(3)={-1.29853,0.580736,0.272342};
Point(4)={-1.13571,0.554571,0.272342};
Point(5)={-0.792725,0.4916,0.272342};
Point(6)={-0.484863,0.430503,0.272342};
Point(7)={-0.229412,0.380106,0.272343};
Point(8)={-0.0213658,0.338825,0.272374};
Point(9)={0.182166,0.296617,0.272146};
Point(10)={0.308832,0.27251,0.27303};
Point(11)={0.36574,0.262206,0.273298};
Point(12)={0.391067,0.256625,0.273682};
Point(13)={0.403053,0.253658,0.273862};
Point(14)={0.410932,0.251707,0.27398};
Point(15)={0.416985,0.2497,0.27406};
Point(16)={0.422054,0.246474,0.274096};
Point(17)={0.427519,0.242996,0.274134};
Point(18)={0.43349,0.239197,0.274176};
Point(19)={0.427515,0.233014,0.27408};
Point(20)={0.41808,0.225635,0.273941};
Point(21)={0.409582,0.222186,0.273819};
Point(22)={0.396835,0.217727,0.273635};
Point(23)={0.373246,0.209727,0.273298};
Point(24)={0.32022,0.193579,0.272711};
Point(25)={0.201697,0.158604,0.272441};
Point(26)={0.0448136,0.112052,0.272278};
Point(27)={-0.230941,0.0297025,0.272331};
Point(28)={-0.547279,-0.064652,0.272342};
Point(29)={-0.953766,-0.186002,0.272342};
Point(30)={-1.61789,-0.382499,0.272342};
Point(31)={-2.62068,-0.659195,0.272342};
Point(32)={-3.78251,-0.931158,0.272342};
Point(33)={-5.14251,-1.09744,0.272342};
Point(34)={-6.42272,-1.19901,0.272342};
Point(35)={-7.74454,-1.22573,0.272342};
Point(36)={-9.27349,-1.22395,0.272342};
Point(37)={-11.0922,-1.22445,0.272342};
Point(38)={-12.8534,-1.22505,0.272342};
Point(39)={-14.2776,-1.22305,0.272342};
Point(40)={-14.9664,-1.22074,0.272342};
Point(41)={-15.5846,-1.19915,0.272342};
Point(42)={-16.3429,-1.14916,0.272342};
Point(43)={-16.8063,-1.10228,0.272342};
Point(44)={-17.1333,-1.05545,0.272342};
Point(45)={-17.4865,-0.981887,0.272342};
Point(46)={-17.7336,-0.904537,0.272342};
Point(47)={-17.852,-0.856089,0.272342};
Point(48)={-17.91,-0.825877,0.272342};
Point(49)={-17.948,-0.79958,0.272342};
Point(50)={-17.9739,-0.776375,0.272342};
Point(51)={-17.9919,-0.754153,0.272342};
Point(52)={-18.005,-0.735759,0.272342};
Point(53)={-18.0195,-0.708138,0.272342};
Point(54)={-18.0233,-0.691911,0.272342};
Point(55)={-18.0173,-0.677519,0.272342};
Point(56)={-18.011,-0.662416,0.272342};
Point(57)={-18.0012,-0.647181,0.272342};
Point(58)={-17.9883,-0.629481,0.272342};
Point(59)={-17.9712,-0.607414,0.272342};
Point(60)={-17.9451,-0.582988,0.272342};
Point(61)={-17.9098,-0.553899,0.272342};
Point(62)={-17.8526,-0.51337,0.272342};
Point(63)={-17.7336,-0.44185,0.272342};
Point(64)={-17.5138,-0.332216,0.272342};
Point(65)={-17.0624,-0.137882,0.272342};
Point(66)={-16.5619,0.0804761,0.272342};
Point(67)={-16.2516,0.268402,0.272342};
Point(68)={-15.9451,0.487979,0.272342};
Point(69)={-15.7811,0.584463,0.272342};
Point(70)={-15.6456,0.640487,0.272342};
Point(71)={-15.4716,0.688214,0.272342};
Point(72)={-15.1592,0.733963,0.272342};
Point(73)={-14.8349,0.753442,0.272342};
Point(74)={-14.2354,0.769157,0.272342};
Point(75)={-12.9053,0.800136,0.272342};
Point(76)={-11.3349,0.804831,0.272342};
Point(77)={-9.7617,0.806311,0.272342};
Point(78)={-8.22982,0.804744,0.272342};
Point(79)={-6.656,0.810412,0.272342};
Point(80)={-5.68102,0.828845,0.272342};
Point(81)={-5.24382,0.842186,0.272342};
Point(82)={-5.04769,0.847178,0.272342};
Point(83)={-4.9597,0.849183,0.272342};
Point(84)={-4.92018,0.850113,0.272342};
Point(85)={-4.8992,0.850588,0.272342};
Point(86)={-4.88341,0.850946,0.272342};
Point(87)={-4.87041,0.864819,0.272342};
Point(88)={-4.84122,0.89556,0.272342};
Point(89)={-4.78676,0.952464,0.272342};
Point(90)={-4.74792,0.993208,0.272342};
Point(91)={-4.71107,1.03196,0.272342};
Point(92)={-4.67099,1.07459,0.272342};
Point(93)={-4.64458,1.10218,0.272342};
Point(94)={-4.61792,1.12984,0.272342};
Point(95)={-4.58084,1.16811,0.272342};
Point(96)={-4.54477,1.20649,0.272342};
Point(97)={-4.51367,1.23922,0.272342};
Point(98)={-4.48381,1.27003,0.272342};
Point(99)={-4.45448,1.30081,0.272342};
Point(100)={-4.42287,1.33375,0.272342};
Point(101)={-4.39,1.3682,0.272342};
Point(102)={-4.35712,1.40303,0.272342};
Point(103)={-4.32433,1.43746,0.272342};
Point(104)={-4.29254,1.47086,0.272342};
Point(105)={-4.25606,1.50876,0.272342};
Point(106)={-4.21472,1.55212,0.272342};
Point(107)={-4.17329,1.59599,0.272342};
Point(108)={-4.13352,1.63741,0.272342};
Point(109)={-4.0899,1.68353,0.272342};
Point(110)={-4.0314,1.74444,0.272342};
Point(111)={-3.99245,1.78471,0.272342};
Point(112)={-3.9553,1.8242,0.272342};
Point(113)={-3.91819,1.8627,0.272342};
Point(114)={-3.8752,1.90762,0.272342};
Point(115)={-3.83256,1.95276,0.272342};
Point(116)={-3.78555,2.00174,0.272342};
Point(117)={-3.75214,2.03709,0.272342};
Point(118)={-3.71593,2.07433,0.272342};
Point(119)={-3.68021,2.11201,0.272342};
Point(120)={-3.64497,2.14858,0.272342};
Point(121)={-3.61207,2.18356,0.272342};
Point(122)={-3.57912,2.21784,0.272342};
Point(123)={-3.53941,2.26024,0.272342};
Point(124)={-3.50377,2.29761,0.272342};
Point(125)={-3.45887,2.34452,0.272342};
Point(126)={-3.41638,2.38884,0.272342};
Point(127)={-3.37656,2.43049,0.272342};
Point(128)={-3.33055,2.47914,0.272342};
Point(129)={-3.28112,2.53071,0.272342};
Point(130)={-3.24876,2.56465,0.272342};
Point(131)={-3.20934,2.60581,0.272342};
Point(132)={-3.16553,2.65221,0.272342};
Point(133)={-3.12681,2.69275,0.272342};
Point(134)={-3.08118,2.74028,0.272342};
Point(135)={-3.02702,2.79722,0.272342};
Point(136)={-2.97475,2.85234,0.272342};
Point(137)={-2.92345,2.90633,0.272342};
Point(138)={-2.86747,2.96449,0.272342};
Point(139)={-2.81318,3.02168,0.272342};
Point(140)={-2.75959,3.0786,0.272342};
Point(141)={-2.71163,3.12868,0.272342};
Point(142)={-2.67168,3.17072,0.272342};
Point(143)={-2.63046,3.21411,0.272342};
Point(144)={-2.58793,3.25859,0.272342};
Point(145)={-2.55052,3.29804,0.272342};
Point(146)={-2.52207,3.32756,0.272342};
Point(147)={-2.4828,3.36864,0.272342};
Point(148)={-2.4343,3.41952,0.272342};
Point(149)={-2.38645,3.46983,0.272342};
Point(150)={-2.33256,3.52633,0.272342};
Point(151)={-2.27413,3.58761,0.272342};
Point(152)={-2.22228,3.64134,0.272342};
Point(153)={-2.178,3.68836,0.272342};
Point(154)={-2.12787,3.74053,0.272342};
Point(155)={-2.07789,3.79274,0.272342};
Point(156)={-2.02957,3.84406,0.272342};
Point(157)={-1.97372,3.90273,0.272342};
Point(158)={-1.92992,3.94846,0.272342};
Point(159)={-1.88092,3.99935,0.272342};
Point(160)={-1.82036,4.06336,0.272342};
Point(161)={-1.76323,4.12278,0.272342};
Point(162)={-1.7114,4.17793,0.272342};
Point(163)={-1.64757,4.24403,0.272342};
Point(164)={-1.60007,4.29422,0.272342};
Point(165)={-1.52382,4.37395,0.272342};
Point(166)={-1.45214,4.44914,0.272342};
Point(167)={-1.39744,4.50571,0.272342};
Point(168)={-1.37216,4.53253,0.272342};
Point(169)={-1.36041,4.54505,0.272342};
Point(170)={-1.32955,4.54505,0.272342};
Point(171)={-1.29567,4.54505,0.272342};
Point(172)={-1.25049,4.54505,0.272342};
Point(173)={-1.15698,4.54505,0.272342};
Point(174)={-1.04935,4.54505,0.272342};
Point(175)={-0.941663,4.54505,0.272342};
Point(176)={-0.835952,4.54505,0.272342};
Point(177)={-0.731249,4.54505,0.272342};
Point(178)={-0.668802,4.54505,0.272342};
Point(179)={-0.601791,4.54505,0.272342};
Point(180)={-0.514994,4.54505,0.272342};
Point(181)={-0.389406,4.54505,0.272342};
Point(182)={-0.314485,4.54505,0.272342};
Point(183)={-0.264262,4.54505,0.272342};
Point(184)={-0.210411,4.54505,0.272342};
Point(185)={-0.216215,4.52553,0.272342};
Point(186)={-0.222808,4.50336,0.272342};
Point(187)={-0.230781,4.47685,0.272342};
Point(188)={-0.245778,4.42983,0.272342};
Point(189)={-0.278844,4.32262,0.272342};
Point(190)={-0.352395,4.07992,0.272342};
Point(191)={-0.521606,3.52739,0.272342};
Point(192)={-0.838727,2.48851,0.272342};
Point(193)={-1.02125,1.89281,0.272342};
Point(194)={-1.14473,1.48527,0.272342};
Point(195)={-1.29254,1.00409,0.272342};
Point(196)={-1.35885,0.787701,0.272342};
Point(197)={-1.38813,0.692385,0.272342};
Point(198)={-1.40116,0.649369,0.272342};
Point(199)={-1.40942,0.622107,0.272342};
Point(200)={-1.41641,0.599016,0.272342};

Line(203)={2,3};
Line(204)={3,4};
Line(205)={4,5};
Line(206)={5,6};
Line(207)={6,7};
Line(208)={7,8};
Line(209)={8,9};
Line(210)={9,10};
Line(211)={10,11};
Line(212)={11,12};
Line(213)={12,13};
Line(214)={13,14};
Line(215)={14,15};
Line(216)={15,16};
Line(217)={16,17};
Line(218)={17,18};
Line(219)={18,19};
Line(220)={19,20};
Line(221)={20,21};
Line(222)={21,22};
Line(223)={22,23};
Line(224)={23,24};
Line(225)={24,25};
Line(226)={25,26};
Line(227)={26,27};
Line(228)={27,28};
Line(229)={28,29};
Line(230)={29,30};
Line(231)={30,31};
Line(232)={31,32};
Line(233)={32,33};
Line(234)={33,34};
Line(235)={34,35};
Line(236)={35,36};
Line(237)={36,37};
Line(238)={37,38};
Line(239)={38,39};
Line(240)={39,40};
Line(241)={40,41};
Line(242)={41,42};
Line(243)={42,43};
Line(244)={43,44};
Line(245)={44,45};
Line(246)={45,46};
Line(247)={46,47};
Line(248)={47,48};
Line(249)={48,49};
Line(250)={49,50};
Line(251)={50,51};
Line(252)={51,52};
Line(253)={52,53};
Line(254)={53,54};
Line(255)={54,55};
Line(256)={55,56};
Line(257)={56,57};
Line(258)={57,58};
Line(259)={58,59};
Line(260)={59,60};
Line(261)={60,61};
Line(262)={61,62};
Line(263)={62,63};
Line(264)={63,64};
Line(265)={64,65};
Line(266)={65,66};
Line(267)={66,67};
Line(268)={67,68};
Line(269)={68,69};
Line(270)={69,70};
Line(271)={70,71};
Line(272)={71,72};
Line(273)={72,73};
Line(274)={73,74};
Line(275)={74,75};
Line(276)={75,76};
Line(277)={76,77};
Line(278)={77,78};
Line(279)={78,79};
Line(280)={79,80};
Line(281)={80,81};
Line(282)={81,82};
Line(283)={82,83};
Line(284)={83,84};
Line(285)={84,85};
Line(286)={85,86};
Line(287)={86,87};
Line(288)={87,88};
Line(289)={88,89};
Line(290)={89,90};
Line(291)={90,91};
Line(292)={91,92};
Line(293)={92,93};
Line(294)={93,94};
Line(295)={94,95};
Line(296)={95,96};
Line(297)={96,97};
Line(298)={97,98};
Line(299)={98,99};
Line(300)={99,100};
Line(301)={100,101};
Line(302)={101,102};
Line(303)={102,103};
Line(304)={103,104};
Line(305)={104,105};
Line(306)={105,106};
Line(307)={106,107};
Line(308)={107,108};
Line(309)={108,109};
Line(310)={109,110};
Line(311)={110,111};
Line(312)={111,112};
Line(313)={112,113};
Line(314)={113,114};
Line(315)={114,115};
Line(316)={115,116};
Line(317)={116,117};
Line(318)={117,118};
Line(319)={118,119};
Line(320)={119,120};
Line(321)={120,121};
Line(322)={121,122};
Line(323)={122,123};
Line(324)={123,124};
Line(325)={124,125};
Line(326)={125,126};
Line(327)={126,127};
Line(328)={127,128};
Line(329)={128,129};
Line(330)={129,130};
Line(331)={130,131};
Line(332)={131,132};
Line(333)={132,133};
Line(334)={133,134};
Line(335)={134,135};
Line(336)={135,136};
Line(337)={136,137};
Line(338)={137,138};
Line(339)={138,139};
Line(340)={139,140};
Line(341)={140,141};
Line(342)={141,142};
Line(343)={142,143};
Line(344)={143,144};
Line(345)={144,145};
Line(346)={145,146};
Line(347)={146,147};
Line(348)={147,148};
Line(349)={148,149};
Line(350)={149,150};
Line(351)={150,151};
Line(352)={151,152};
Line(353)={152,153};
Line(354)={153,154};
Line(355)={154,155};
Line(356)={155,156};
Line(357)={156,157};
Line(358)={157,158};
Line(359)={158,159};
Line(360)={159,160};
Line(361)={160,161};
Line(362)={161,162};
Line(363)={162,163};
Line(364)={163,164};
Line(365)={164,165};
Line(366)={165,166};
Line(367)={166,167};
Line(368)={167,168};
Line(369)={168,169};
Line(370)={169,170};
Line(371)={170,171};
Line(372)={171,172};
Line(373)={172,173};
Line(374)={173,174};
Line(375)={174,175};
Line(376)={175,176};
Line(377)={176,177};
Line(378)={177,178};
Line(379)={178,179};
Line(380)={179,180};
Line(381)={180,181};
Line(382)={181,182};
Line(383)={182,183};
Line(384)={183,184};
Line(385)={184,185};
Line(386)={185,186};
Line(387)={186,187};
Line(388)={187,188};
Line(389)={188,189};
Line(390)={189,190};
Line(391)={190,191};
Line(392)={191,192};
Line(393)={192,193};
Line(394)={193,194};
Line(395)={194,195};
Line(396)={195,196};
Line(397)={196,197};
Line(398)={197,198};
Line(399)={198,199};
Line(400)={199,200};
Line(401)={200,2};



// Mid plane
x_dim[] = {-20.0,1.0};
y_dim[] = {-2.0,5.0};
zv=0.272334214543;

mp1 = newp; Point(mp1) = {x_dim[0], y_dim[0], zv};
mp2 = newp; Point(mp2) = {x_dim[1], y_dim[0], zv};
mp3 = newp; Point(mp3) = {x_dim[1], y_dim[1], zv};
mp4 = newp; Point(mp4) = {x_dim[0], y_dim[1], zv};

ml1 = newl; Line(ml1) = {mp1,mp2};
ml2 = newl; Line(ml2) = {mp2,mp3};
ml3 = newl; Line(ml3) = {mp3,mp4};
ml4 = newl; Line(ml4) = {mp4,mp1};

mll  = newll; Line Loop(mll) = {ml1,  ml2, ml3, ml4};
mplane  = news; Plane Surface(mplane) = {mll}; // Bottom unstructured region

mplane_ribs[] = {ml1,  ml2, ml3, ml4};
mid_plane[] = {mplane};

internal_ribs[]={203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,277,278,279,280,281,282,283,284,285,286,287,288,289,290,291,292,293,294,295,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,314,315,316,317,318,319,320,321,322,323,324,325,326,327,328,329,330,331,332,333,334,335,336,337,338,339,340,341,342,343,344,345,346,347,348,349,350,351,352,353,354,355,356,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385,386,387,388,389,390,391,392,393,394,395,396,397,398,399,400,401};

Physical Line("InternalRibs") = {internal_ribs[]};
Physical Line("ExtenalRibs") = {mplane_ribs[]};
imfalconoll  = newll; Line Loop(imfalconoll) = {internal_ribs[]};
//omfalconill  = newll; Line Loop(omfalconill) = {mplane_ribs[]};

falconplane  = news; Surface(falconplane) = {imfalconoll}; // Bottom unstructured region


