      SUBROUTINE COEFFT
      include 'PARAM.include'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/ARRAYS/S(NBSZR,NBSZR),Y(9135),Z(765),XX(NAIGAIO1)
      DO 10 I=1,9135
   10 Y(I)=0.0D0
      DO 20 I=1,765
   20 Z(I)=0.0D0
C     LOAD NON-ZERO Y COEFFICIENTS
      Y(7039)=    64.D0 
      Y(7040)=    64.D0 
      Y(7049)=   -64.D0 
      Y(7032)=  -128.D0 
      Y(7041)=   -64.D0 
      Y(7033)=  -128.D0 
      Y(7042)=   128.D0 
      Y(7025)=    64.D0 
      Y(7034)=   128.D0 
      Y(7026)=    64.D0 
      Y(7035)=   -64.D0 
      Y(7027)=   -64.D0 
      Y(6904)=   -96.D0 
      Y(6913)=    32.D0 
      Y(6896)=  -192.D0 
      Y(6905)=   192.D0 
      Y(6906)=   288.D0 
      Y(6915)=   -96.D0 
      Y(6889)=   192.D0 
      Y(6907)=  -192.D0 
      Y(6890)=    96.D0 
      Y(6899)=  -288.D0 
      Y(6891)=  -192.D0 
      Y(6900)=   192.D0 
      Y(6892)=   -32.D0 
      Y(6901)=    96.D0 
      Y(2854)=   -16.D0 
      Y(2863)=    16.D0 
      Y(2847)=    32.D0 
      Y(2856)=   -16.D0 
      Y(2865)=   -16.D0 
      Y(2840)=   -16.D0 
      Y(2849)=   -16.D0 
      Y(2858)=    32.D0 
      Y(2842)=    16.D0 
      Y(2851)=   -16.D0 
      Y(2710)=    48.D0 
      Y(2719)=   -48.D0 
      Y(2711)=    48.D0 
      Y(2720)=   -96.D0 
      Y(2729)=    48.D0 
      Y(2703)=   -48.D0 
      Y(2712)=   -48.D0 
      Y(2721)=    96.D0 
      Y(2704)=   -48.D0 
      Y(2713)=    48.D0 
      Y(2722)=    48.D0 
      Y(2731)=   -48.D0 
      Y(2705)=    96.D0 
      Y(2714)=   -48.D0 
      Y(2723)=   -48.D0 
      Y(2706)=    48.D0 
      Y(2715)=   -96.D0 
      Y(2724)=    48.D0 
      Y(2707)=   -48.D0 
      Y(2716)=    48.D0 
      Y(5329)=    64.D0 
      Y(5322)=  -128.D0 
      Y(5340)=   -64.D0 
      Y(5315)=    64.D0 
      Y(5333)=   128.D0 
      Y(5326)=   -64.D0 
      Y(5185)=   -96.D0 
      Y(5194)=    32.D0 
      Y(5186)=   -96.D0 
      Y(5195)=    64.D0 
      Y(5204)=    32.D0 
      Y(5178)=    96.D0 
      Y(5187)=    32.D0 
      Y(5196)=    64.D0 
      Y(5179)=    96.D0 
      Y(5188)=   -32.D0 
      Y(5197)=    32.D0 
      Y(5206)=   -96.D0 
      Y(5180)=   -64.D0 
      Y(5189)=   -32.D0 
      Y(5198)=   -96.D0 
      Y(5181)=   -32.D0 
      Y(5190)=   -64.D0 
      Y(5199)=    96.D0 
      Y(5182)=   -32.D0 
      Y(5191)=    96.D0 
      Y(4375)=  -144.D0 
      Y(4384)=    96.D0 
      Y(4393)=   -16.D0 
      Y(4368)=   144.D0 
      Y(4386)=   -48.D0 
      Y(4395)=    96.D0 
      Y(4370)=   -96.D0 
      Y(4379)=    48.D0 
      Y(4397)=  -144.D0 
      Y(4372)=    16.D0 
      Y(4381)=   -96.D0 
      Y(4390)=   144.D0 
      Y(1900)=   144.D0 
      Y(1909)=  -144.D0 
      Y(1893)=  -144.D0 
      Y(1920)=   144.D0 
      Y(1895)=   144.D0 
      Y(1922)=  -144.D0 
      Y(1906)=  -144.D0 
      Y(1915)=   144.D0 
      Y( 955)=   -16.D0 
      Y( 964)=    32.D0 
      Y( 973)=   -16.D0 
      Y( 948)=    16.D0 
      Y( 966)=   -48.D0 
      Y( 975)=    32.D0 
      Y( 959)=    48.D0 
      Y( 977)=   -16.D0 
      Y( 952)=    16.D0 
      Y( 961)=   -32.D0 
      Y( 970)=    16.D0 
      Y(8155)=    64.D0 
      Y(8156)=   -64.D0 
      Y(8165)=   -64.D0 
      Y(8148)=   -64.D0 
      Y(8157)=    64.D0 
      Y(8149)=    64.D0 
      Y(8158)=    64.D0 
      Y(8150)=   -64.D0 
      Y(8020)=   -96.D0 
      Y(8029)=    32.D0 
      Y(8021)=   128.D0 
      Y(8013)=    96.D0 
      Y(8031)=   -96.D0 
      Y(8014)=  -128.D0 
      Y(8015)=   -32.D0 
      Y(8024)=    96.D0 
      Y(7084)=   -64.D0 
      Y(7076)=  -128.D0 
      Y(7085)=    64.D0 
      Y(7086)=   128.D0 
      Y(7069)=   128.D0 
      Y(7070)=    64.D0 
      Y(3214)=    16.D0 
      Y(3205)=   -16.D0 
      Y(7071)=   -64.D0 
      Y(7079)=  -128.D0 
      Y(3206)=    16.D0 
      Y(3215)=   -16.D0 
      Y(3198)=    16.D0 
      Y(3216)=   -16.D0 
      Y(3199)=   -16.D0 
      Y(3217)=    16.D0 
      Y(3200)=   -16.D0 
      Y(3209)=    16.D0 
      Y(3201)=    16.D0 
      Y(3210)=   -16.D0 
      Y(7579)=    64.D0 
      Y(7580)=   -64.D0 
      Y(7572)=  -128.D0 
      Y(7573)=   128.D0 
      Y(7565)=    64.D0 
      Y(7566)=   -64.D0 
      Y(5680)=    64.D0 
      Y(5681)=   -64.D0 
      Y(5673)=   -64.D0 
      Y(5691)=   -64.D0 
      Y(5674)=    64.D0 
      Y(5692)=    64.D0 
      Y(5684)=    64.D0 
      Y(5685)=   -64.D0 
      Y(7435)=   -96.D0 
      Y(7444)=    32.D0 
      Y(7436)=   -96.D0 
      Y(7445)=   160.D0 
      Y(7428)=    96.D0 
      Y(7437)=   128.D0 
      Y(7446)=   -96.D0 
      Y(7429)=    96.D0 
      Y(7438)=  -128.D0 
      Y(7447)=   -96.D0 
      Y(7430)=  -160.D0 
      Y(7439)=    96.D0 
      Y(7431)=   -32.D0 
      Y(7440)=    96.D0 
      Y(5545)=   -96.D0 
      Y(5554)=    32.D0 
      Y(5546)=    32.D0 
      Y(5555)=    32.D0 
      Y(5538)=    96.D0 
      Y(5556)=    32.D0 
      Y(5539)=   -32.D0 
      Y(5557)=   -96.D0 
      Y(5540)=   -32.D0 
      Y(5549)=   -32.D0 
      Y(5541)=   -32.D0 
      Y(5550)=    96.D0 
      Y(3079)=   -48.D0 
      Y(3070)=    48.D0 
      Y(3071)=   -48.D0 
      Y(3080)=    48.D0 
      Y(3063)=   -48.D0 
      Y(3081)=    48.D0 
      Y(3064)=    48.D0 
      Y(3082)=   -48.D0 
      Y(3065)=    48.D0 
      Y(3074)=   -48.D0 
      Y(3066)=   -48.D0 
      Y(3075)=    48.D0 
      Y(8200)=   -64.D0 
      Y(8201)=    64.D0 
      Y(8193)=    64.D0 
      Y(8194)=   -64.D0 
      Y(7615)=   -64.D0 
      Y(7616)=   -64.D0 
      Y(7625)=    64.D0 
      Y(7608)=    64.D0 
      Y(7617)=    64.D0 
      Y(7609)=    64.D0 
      Y(7618)=   -64.D0 
      Y(7610)=   -64.D0 
      Y(3250)=    16.D0 
      Y(3259)=   -16.D0 
      Y(3243)=   -16.D0 
      Y(3261)=    16.D0 
      Y(3245)=    16.D0 
      Y(3254)=   -16.D0 
      Y(5725)=   -64.D0 
      Y(5718)=    64.D0 
      Y(5736)=    64.D0 
      Y(5729)=   -64.D0 
C     LOAD NON-ZERO Z COEFFICIENTS
      Z(341)=    -1.D0 
      Z(343)=     3.D0 
      Z(345)=    -3.D0 
      Z(347)=     1.D0 
      Z(664)=    -1.D0 
      Z(665)=     5.D0 
      Z(666)=   -10.D0 
      Z(667)=    10.D0 
      Z(668)=    -5.D0 
      Z(669)=     1.D0 
      Z(154)=    -1.D0 
      Z(156)=     5.D0 
      Z(158)=   -10.D0 
      Z(160)=    10.D0 
      Z(162)=    -5.D0 
      Z(164)=     1.D0 
      Z(222)=    -1.D0 
      Z(223)=     1.D0 
      Z(224)=     4.D0 
      Z(225)=    -4.D0 
      Z(226)=    -6.D0 
      Z(227)=     6.D0 
      Z(228)=     4.D0 
      Z(229)=    -4.D0 
      Z(230)=    -1.D0 
      Z(231)=     1.D0 
      Z(307)=    -1.D0 
      Z(308)=     2.D0 
      Z(309)=     2.D0 
      Z(310)=    -6.D0 
      Z(312)=     6.D0 
      Z(313)=    -2.D0 
      Z(314)=    -2.D0 
      Z(315)=     1.D0 
      Z(409)=    -1.D0 
      Z(410)=     3.D0 
      Z(411)=    -1.D0 
      Z(412)=    -5.D0 
      Z(413)=     5.D0 
      Z(414)=     1.D0 
      Z(415)=    -3.D0 
      Z(416)=     1.D0 
      Z(528)=    -1.D0 
      Z(529)=     4.D0 
      Z(530)=    -5.D0 
      Z(532)=     5.D0 
      Z(533)=    -4.D0 
      Z(534)=     1.D0 
      Z(562)=    -1.D0 
      Z(563)=     2.D0 
      Z(565)=    -2.D0 
      Z(566)=     1.D0 
      Z(732)=    -1.D0 
      Z(733)=     1.D0 
      Z(545)=     1.D0 
      Z(546)=    -3.D0 
      Z(547)=     2.D0 
      Z(548)=     2.D0 
      Z(549)=    -3.D0 
      Z(550)=     1.D0 
      Z(579)=     1.D0 
      Z(580)=    -1.D0 
      Z(581)=    -1.D0 
      Z(582)=     1.D0 
      Z(596)=    -1.D0 
      Z(598)=     1.D0 
      Z(443)=    -1.D0 
      Z(444)=     1.D0 
      Z(445)=     2.D0 
      Z(446)=    -2.D0 
      Z(447)=    -1.D0 
      Z(448)=     1.D0 
      Z(698)=    -1.D0 
      Z(699)=     3.D0 
      Z(700)=    -3.D0 
      Z(701)=     1.D0 
      Z(324)=     1.D0 
      Z(325)=    -1.D0 
      Z(326)=    -3.D0 
      Z(327)=     3.D0 
      Z(328)=     3.D0 
      Z(329)=    -3.D0 
      Z(330)=    -1.D0 
      Z(331)=     1.D0 
      Z(460)=     1.D0 
      Z(462)=    -2.D0 
      Z(464)=     1.D0 
      RETURN
      END