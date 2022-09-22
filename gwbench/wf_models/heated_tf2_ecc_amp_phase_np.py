# Copyright (C) 2020  Ssohrab Borhanian
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


import numpy as np

import gwbench.basic_relations as brs
from gwbench.basic_constants import time_fac, strain_fac

cos = np.cos
sin = np.sin
log = np.log
ustep = np.heaviside
PI = np.pi

wf_symbs_string = 'f Mc eta chi1z chi2z DL tc phic iota Heff5 Heff8 e0'

#------from Anruadha--------
# defining constants
GammaE = 0.577215664901532

def hfpc(f, Mc, eta, chi1z, chi2z, DL, tc, phic, iota, Heff5, Heff8, e0):
    '''
    Mc ... in solar mass
    DL ... in mega parsec
    '''
    # convert to sec
    Mc = Mc * time_fac
    DL = DL * time_fac/strain_fac

    # get sym and asym chi combinations
    chi_s = brs.chi_s(chi1z,chi2z)
    chi_a = brs.chi_a(chi1z,chi2z)

    '''
    Mc is in sec, e.g., Mc = 10*MTSUN_SI (for 10 solar mass)
    DL is in sec, e.g., DL = 100*1e6*PC_SI/C_SI (for 100 Mpc)
    '''

    M = Mc/eta**(3./5.)
    delta = (1.-4.*eta)**0.5
    ci = cos(iota)
    si = sin(iota)
    v  = (PI*M*f)**(1./3.)
    v1 = (PI*M*2*f/1)**(1./3.)
    v2 = (PI*M*2*f/2)**(1./3.)
    v3 = (PI*M*2*f/3)**(1./3.)
    v4 = (PI*M*2*f/4)**(1./3.)
    v5 = (PI*M*2*f/5)**(1./3.)
    v6 = (PI*M*2*f/6)**(1./3.)
    v7 = (PI*M*2*f/7)**(1./3.)
    v8 = (PI*M*2*f/8)**(1./3.)
    flso = brs.f_isco(M)
    vlso = (PI*M*flso)**(1./3.)
    f0 = f[0]
    v0 = (PI*M*f0)**(1./3.)
    vByv0 = v/v0
    A =((5./24.)**0.5/PI**(2./3.))*(Mc**(5./6.)/DL)  #Amplitude for the eccentric case
   # A =((5./24.)**0.5/PI**(2./3.))*(Mc**(5./6.)/DL) Amplitude for the circular case

  # Sixth order eccentricity contributions in leading order GW amplitude: Plus polarization. We set the azimuthal angle beta as zero. 
    hp1 = ((3./2.)*(1+ci**2.) + si**2.)*vByv0**(-19./6.)*e0 + (((-30299./3648.)*(1+ci**2.) - (9517./1824.)*si**2.)*vByv0**(-19./2.) + ((3323./1216.)*(1+ci**2.) + (3323./1824.)*si**2.)*vByv0**(-19./6.))*e0**3. + (((717415013./13307904.)*(1+ci**2.) + (220389695./6653952.)*si**2.)*vByv0**(-95./6.) + ((-100683577./2217984.)*(1+ci**2.) - (31624991./1108992.)*si**2.)*vByv0**(-19./2.) + ((15994231./4435968.)*(1+ci**2.) + (15994231./6653952.)*si**2.)*vByv0**(-19./6.))*e0**5.
    hp2 = (-2.)*(1+ci**2.) + ((277./24.)*(1+ci**2.) + si**2.)*vByv0**(-19./3.)*e0**2. + (((-3254599./43776.)*(1+ci**2.) - (3305./456.)*si**2.)*vByv0**(-38./3.) + ((920471./21888.)*(1+ci**2.) + (3323./912)*si**2.)*vByv0**(-19./3.))*e0**4. + (((103729904239./199618560.)*(1+ci**2.) + (29064841./554496.)*si**2.)*vByv0**(-19.) + ((-10815032477./19961856.)*(1+ci**2.) - (10982515./207936.)*si**2.)*vByv0**(-38./3.) + ((468070445./4990464.)*(1+ci**2.) + (1689785./207936.)*si**2.)*vByv0**(-19./3.))*e0**6.
    hp3 = (-9./2.)*(1+ci**2.)*vByv0**(-19./6.)*e0 + (((40863./1216.)*(1+ci**2.) + (9./8.)*si**2.)*vByv0**(-19./2.) - (9969./1216.)*(1+ci**2.)*vByv0**(-19./6.))*e0**3. + (((-1810486747./7393280.)*(1+ci**2.) - (50883./4864.)*si**2.)*vByv0**(-95./6.) + ((135787749./739328.)*(1+ci**2.) + (29907./4864.)*si**2.)*vByv0**(-19./2.) - (15994231./1478656.)*(1+ci**2.)*vByv0**(-19./6.))*e0**5.
    hp4 = (-8.)*(1+ci**2.)*vByv0**(-19./3.)*e0**2. + (((1431./19.)*(1+ci**2.) + (4./3.)*si**2.)*vByv0**(-38./3.) - (3323./114.)*(1+ci**2.)*vByv0**(-19./3.))*e0**4. + (((-71474367./115520.)*(1+ci**2.) - (51793./3420.)*si**2.)*vByv0**(-19.) + ((1585071./2888.)*(1+ci**2.) + (3323./342.)*si**2.)*vByv0**(-38./3.) - (1689785./25992.)*(1+ci**2.)*vByv0**(-19./3.))*e0**6.
    hp5 = (-625./48.)*(1+ci**2.)*vByv0**(-19./2.)*e0**3. + (((13023125./87552.)*(1+ci**2.) + (625./384.)*si**2.)*vByv0**(-95./6.) - (2076875./29184.)*(1+ci**2.)*vByv0**(-19./2.))*e0**5.
    hp6 = (-81./4.)*(1+ci**2.)*vByv0**(-38./3.)*e0**4. + (((1656963./6080.)*(1+ci**2.) + (81./40.)*si**2.)*vByv0**(-19.) - (89721./608.)*(1+ci**2.)*vByv0**(-38./3.))*e0**6.
    hp7 = (-117649./3840.)*(1+ci**2.)*vByv0**(-95./6.)*e0**5.
    hp8 = (-2048./45.)*(1+ci**2.)*vByv0**(-19.)*e0**6.
	
  # Sixth order eccentricity contributions in leading order GW amplitude: Cross polarization. We set the azimuthal angle beta as zero. 
    hc1 =  3.j*ci*vByv0**(-19./6.)*e0. + ((-31363./1824.)*vByv0**(-19./2.) + (3323./608.)*vByv0**(-19./6.))*1j*ci*e0**3. + ((749695861./6653952.)*vByv0**(-95./6.) - (104219249./1108992.)*vByv0**(-19./2.) + (15994231./2217984.)*vByv0**(-19./6.))*1j*ci*e0**5.
    hc2 = -4.j*ci + (277.j*ci/12.)*vByv0**(-19./3.)*e0**2. + ((-3265543./21888.)*vByv0**(-38./3.) + (920471./10944.)*vByv0**(-19./3.))*1j*ci*e0**4. + ((104238504751./99809280.)*vByv0**(-19.) - (10851399389./9980928.)*vByv0**(-38./3.) + (468070445./2495232.)*vByv0**(-19./3.))*1j*ci*e0**6.
    hc3 = (-9.j*ci)*vByv0**(-19./6.)*e0 + ((40863./608.)*vByv0**(-19./2.) - (9969./608.)*vByv0**(-19./6.))*1j*ci*e0**3. + (-(1812254203./3696640.)*vByv0**(-95./6.) + (135787749./369664.)*vByv0**(-19./2.) - (15994231./739328.)*vByv0**(-19./6.))*1j*ci*e0**5.
    hc4 = (-16.j*ci)*vByv0**(-19./3.)*e0**2. + ((2862./19.)*vByv0**(-38./3.)-(3323./57.)*vByv0**(-19./3.))*1j*ci*e0**4. + (-(643523447./519840.)*vByv0**(-19.) + (1585071./1444.)*vByv0**(-38./3.) - (1689785./12996.)*vByv0**(-19./3.))*1j*ci*e0**6.
    hc5 = (-625.j*ci/24.)*vByv0**(-19./2.)*e0**3. + ((13023125./43776.)*vByv0**(-95./6.) - (2076875./14592.)*vByv0**(-19./2.))*1j*ci*e0**5.
    hc6 = (-81.j*ci/2.)*vByv0**(-38./3.)*e0**4. + ((1656963./3040.)*vByv0**(-19.) - (89721./304.)*vByv0**(-38./3.))*1j*ci*e0**6.
    hc7 = (-117649.j*ci/1920.)*vByv0**(-95./6.)*e0**5.
    hc8 = (-4096.j*ci/45.)*vByv0**(-19.)*e0**6.


  # 3.5PN phasing (point particle limit)
    p0 = 1. 

    p1 = 0

    p2 = (3715./756. + (55.*eta)/9.)

    p3 = (-16.*PI + (113.*delta*chi_a)/3. + (113./3. - (76.*eta)/3.)*chi_s)

    p4 = (15293365./508032. + (27145.*eta)/504.+ (3085.*eta**2)/72. + (-405./8. + 200.*eta)*chi_a**2 - (405.*delta*chi_a*chi_s)/4. + (-405./8. + (5.*eta)/2.)*chi_s**2)

    gamma = (732985./2268. - 24260.*eta/81. - 340.*eta**2/9.)*chi_s + (732985./2268. + 140.*eta/9.)*delta*chi_a

    p5 = (38645.*PI/756. - 65.*PI*eta/9. - gamma)

    p5L = (38645.*PI/756. - 65.*PI*eta/9. - gamma)*3*log(v)  # It should just be log(v) and not log(v/vlso)

    p6 = (11583231236531./4694215680. - 640./3.*PI**2 - 6848./21.*GammaE + eta*(-15737765635./3048192. + 2255./12.*PI**2) + eta*eta*76055./1728. - eta*eta*eta*127825./1296. \
         - (6848./21.)*log(4.) + PI*(2270.*delta*chi_a/3. + (2270./3. - 520.*eta)*chi_s) + (75515./144. - 8225.*eta/18.)*delta*chi_a*chi_s \
         + (75515./288. - 263245.*eta/252. - 480.*eta**2)*chi_a**2 + (75515./288. - 232415.*eta/504. + 1255.*eta**2/9.)*chi_s**2)

    p6L = -(6848./21.)*log(v)

    p7 = (((77096675.*PI)/254016. + (378515.*PI*eta)/1512.- (74045.*PI*eta**2)/756. + (-25150083775./3048192. + (10566655595.*eta)/762048. - (1042165.*eta**2)/3024. + (5345.*eta**3)/36.
         + (14585./8. - 7270.*eta + 80.*eta**2)*chi_a**2)*chi_s + (14585./24. - (475.*eta)/6. + (100.*eta**2)/3.)*chi_s**3 + delta*((-25150083775./3048192.
         + (26804935.*eta)/6048. - (1985.*eta**2)/48.)*chi_a + (14585./24. - 2380.*eta)*chi_a**3 + (14585./8. - (215.*eta)/2.)*chi_a*chi_s**2)))

    # 3PN phasing for point particles on eccentric orbits
    p0ecc = (-(2355./1462.)*vByv0**(-19./3.)*e0**2.) + ((-(2608555./444448.)*vByv0**(-19./3.) + (5222765./998944.)*vByv0**(-38./3.))*e0**4.) + ((-(1326481225./101334144.)*vByv0**(-19./3.) \
	    + (17355248095./455518464.)*vByv0**(-38./3.) - (75356125./3326976.)*vByv0**(-19.))*e0**6.)

    p1ecc = 0

    p2ecc = (((-(2045665./348096.) + (-(128365.*eta)/12432.))*vByv0**(-19./3.) + (-(2223905./491232.) + ((154645.*eta)/17544.))*vByv0**(-25./3.))*e0**2.) + (((-(6797744795./317463552.) - (426556895.*eta/11337984.))*vByv0**(-19./3.) \
            + (-(14275935425./416003328.) + (209699405.*eta/4000032.))*vByv0**(-25./3.) + (198510270125./10484877312. + (1222893635.*eta)/28804608.)*vByv0**(-38./3.) + (14796093245./503467776. - (1028884705.*eta/17980992.))*vByv0**(-44./3.))*e0**4.) \
            + (((-(3456734032025./72381689856.) - (216909251525.*eta/2585060352.))*vByv0**(-19./3.) + (-(2441897241139735./21246121967616.) + (9479155594325.*eta)/58368466944.)*vByv0**(-25./3.) + (659649627625375./4781104054272. \
	    + (4063675549105.*eta)/13134901248.)*vByv0**(-38./3.) + (1968906345873305./5969113952256. - (8999675405695.*eta/16398664704.))*vByv0**(-44./3.) + (-(144936872901./1691582464.) - (7378552295.*eta/32530432.))*vByv0**(-19.) \
	    + (-(213483902125./1117863936.) + (14845156625.*eta)/39923712.)*vByv0**(-21.))*e0**6.)

    p3ecc = (((65561./4080.)*vByv0**(-19./3.) + (-(295945./35088.))*vByv0**(-28./3.))*PI*e0**2.) + (((217859203./3720960.)*vByv0**(-19./3.) - (3048212305./64000512.)*vByv0**(-28./3.) - (6211173025./102085632.)*vByv0**(-38./3.) \
            + (1968982405./35961984.)*vByv0**(-47./3.))*PI*e0**4.) + (((22156798877./169675776.)*vByv0**(-19./3.) - (126468066221755./846342770688.)*vByv0**(-28./3.) - (20639727962075./46551048192.)*vByv0**(-38./3.) + (33366234820475./65594658816.)*vByv0**(-47./3.) \
	    + (30628811474315./97254162432.)*vByv0**(-19.) - (28409259125./79847424.)*vByv0**(-22.))*PI*e0**6.)

    p4ecc = (((-(111064865./14141952.) - (165068815.*eta)/4124736. - (10688155.*eta**2.)/294624.)*vByv0**(-19./3.) + (-(5795368945./350880768.) + (4917245.*eta)/1566432. + (25287905.*eta**2)/447552.)*vByv0**(-25./3.) + (936702035./1485485568. + (3062285.*eta)/260064. \
            - (14251675.*eta**2.)/631584.)*vByv0**(-31./3.))*e0**2.) + (((-(369068546395./12897460224.) - (548523672245.*eta)/3761759232. - (35516739065.*eta**2.)/268697088.)*vByv0**(-19./3.) + (-(37202269351825./297145884672.) - (2132955527705.*eta)/74286471168. \
	    + (34290527545.*eta**2.)/102041856.)*vByv0**(-25./3.) + (-(94372278903235./7251965779968.) + (126823556396665.*eta/733829870592.) - (20940952805.*eta**2./93768192.))*vByv0**(-31./3.) + (418677831611033./34573325230080. + (2163514670909.*eta/12862100160.) \
	    + (203366083643.*eta**2./1130734080.))*vByv0**(-38./3.) + (562379595264125./5284378165248. + (2965713234395.*eta/94363895808.) - (240910046095.*eta**2./518482944.))*vByv0**(-44./3.) + (3654447011975./98224939008. - (4300262795285.*eta/18124839936.) \
	    + (392328884035.*eta**2./1294631424.))*vByv0**(-50./3.))*e0**4.) + ((((-187675742904025./2940620931072.) - (278930807554775.*eta/857681104896.) - (18060683996675.*eta**2./61262936064.))*vByv0**(-19./3.) + ((-6363444229039638215./15175834621968384.) \
	    - (39088433492776445.*eta/270997046820864.) + (1550053258427425.*eta**2/1488994762752.))*vByv0**(-25./3.) + ((-387035983120116605285./5846592827536441344.) + (1095104635088909345.*eta/1338505683959808.) - (185468261986684025.*eta**2./191215097708544.))*vByv0**(-31./3.) 
	    + ((1391266434443462659./15765436304916480.) + (7189359251430607.*eta/5865117672960.) + (675785495945689.*eta**2./515614740480.))*vByv0**(-38./3.) + ((74835480932061169625./62651587527180288.) + (14868442349448515.*eta/21514968244224.) - (2107245064767505*eta**2./472856444928.))*vByv0**(-44./3.) \
	    + ((43949506831840859555./63177102070677504.) - (1344731894414361455.*eta/376054178992128) + (7946157848161165.*eta**2./2066231752704.))*vByv0**(-50./3.) + ((-984783138418096685./40879050017734656.) - (258954290041765.*eta/271268315136.) - (173415564792655.*eta**2./148551696384.))*vByv0**(-19.) \
	    + ((-136868720309511./189457235968.) - (17969188685519.*eta/35523231744.) + (1453574802115.*eta**2./390365184.))*vByv0**(-21.) + ((-26945014260125./52819070976.) + (17350371000625.*eta/6707183616.) - (357715525375.*eta**2./119771136.))*vByv0**(-23.))*e0**6.)

    p5ecc = (((3873451./100548. + (15803101.*eta)/229824.)*vByv0**(-19./3.) + (185734313./4112640. - (12915517.*eta)/146880.)*vByv0**(-25./3.) + (-(771215705 /25062912) - (48393605.*eta)/895104.)*vByv0**(-28./3.) + (-(7063901./520128.) + (149064749.*eta)/2210544.)*vByv0**(-34./3.))*PI*e0**2.) \
	    + ((((12871477673./91699776.) + (52513704623.*eta/209599488.))*vByv0**(-19./3.) + ((238457223541./696563712.) - (17513506613.*eta/33488640.))*vByv0**(-25./3.) + ((-7943466528545./45714751488.) - (498450665645.*eta/1632669696.))*vByv0**(-28./3.) + ((-2408172473789./6790791168.) \
	    + (992200223893.*eta/1697697792.))*vByv0**(-34./3.) + ((6135798516097./2049044594688.) - (15387742160333.*eta/39404703744.))*vByv0**(-38./3.) + ((-17596253179825./51451158528.) + (1223601085925.*eta/1837541376.))*vByv0**(-44./3.) + ((5756797833625./29035044864.) + (461030900395.*eta/1036965888.))*vByv0**(-47./3.) \
	    + ((14896370333./61544448.) - (351697861441.*eta/476969472.))*vByv0**(-53./3.))*PI*e0**4.) + ((((6545299398035./20907548928.) + (26703843023285.*eta/47788683264.))*vByv0**(-19./3.) + ((203940414046321231./177874509496320.) - (158334501890329.*eta/97733246976.))*vByv0**(-25./3.) \
	    + ((-329568530812135595./604531873677312.) - (20680348179051695.*eta/21590424059904.))*vByv0**(-28./3.) + ((-279594780479556044255./145760537338970112.) + (48634782568328640205.*eta/19621610795630592.))*vByv0**(-34./3.) + ((20389258468990331./934364335177728.) - (51133467198786559.*eta/17968544907264.))*vByv0**(-38./3.) \
	    + ((-2341521777112236925./610004935507968.) + (10702863543278075.*eta/1675837734912.))*vByv0**(-44./3.) + ((1268205689374626875./688478983815168.) + (7812596619965525.*eta/1891425779712.))*vByv0**(-47./3.) + ((616055512637722733./132238832173056.) - (292997755491718561.*eta/33059708043264.))*vByv0**(-53./3.) \
	    + ((-8404454631967383095./32413385944203264.) + (1664283962654437115.*eta/623334345080832.))*vByv0**(-19.) + ((86771422906734395./32677398577152.) - (6033875860440055.*eta/1167049949184.))*vByv0**(-21.) + ((-1401056438043./1040973824.) - (2781714215215.*eta/780730368.))*vByv0**(-22.) + ((-34512939466525./13414367232.) \
	    + (22598442827675.*eta/3353591808.))*vByv0**(-24.))*PI*e0**6.)

    p5Lecc = 0

    p6ecc = (((59648637301056877./112661176320000. - (21508213.*PI**2.)/276480. - (734341.*GammaE)/16800. - (409265200567.*eta)/585252864. + (103115.*PI**2.*eta)/6144. - (4726688461.*eta**2.)/34836480. \
	    - (69237581.*eta**3.)/746496. - ((9663919.*log(2.))/50400.) + (4602177.*log(3.))/44800.)*vByv0**(-19./3.) + (-(314646762545./14255087616.) - (1733730575525.*eta)/24946403328. \
	    + (11585856665.*eta**2.)/98993664. + (2105566535.*eta**3.)/10606464.)*vByv0**(-25./3.) + ((24716497.*PI**2.)/293760.)*vByv0**(-28./3.) + (2440991806915./1061063442432. \
            + (1781120054275.*eta)/37895122944. - (1029307085.*eta**2.)/150377472. - (2330466575.*eta**3.)/16111872.)*vByv0**(-31./3.) + (-(4165508390854487./16471063977984.) \
            - (96423905.*PI**2.)/5052672. + (2603845.*GammaE)/61404. - (1437364085977.*eta)/53477480448. + (3121945.*PI**2.*eta)/561408. + (4499991305.*eta**2.)/636636672. + (2425890995.*eta**3.)/68211072. \
            + (1898287.*log(2.))/184212. + (12246471.*log(3.))/163744.)*vByv0**(-37./3.))*e0**2.) + (((198212421751412002271./102746992803840000. - (71471791799.*PI**2./252149760.) - (2440215143.*GammaE/15321600.) \
	    - (1359988261484141.*eta/533750611968.) + (342651145.*PI**2.*eta/5603328.) - (15706785755903.*eta**2./31770869760.) - (230076481663.*eta**3./680804352.) - (32113202837.*log(2.)/45964800.) + (5097678057.*log(3.)/13619200.))*vByv0**(-19./3.) \
	    + ((-2019815083727825./12072022769664.) - (3152945060595815.*eta/5281509961728.) + (34531134933245.*eta**2./65203826688.) + (2855158909615.*eta**3./2418273792.))*vByv0**(-25./3.) + (254578148953*PI**2./535818240.)*vByv0**(-28./3.) \
	    + ((-4180788732081485155./88059777214316544.) + (245803713210749465.*eta/449284577624064.) + (522467227468915.*eta**2./1782875308032.) - (147245442666235.*eta**3./102858190848.))*vByv0**(-31./3.) + ((102453749612934666311./19868699733442560.) \
	    - (598067688595.*PI**2./4608036864.) - (36290762107.*GammaE/56000448.) + (6738669506224179365.*eta/2219101528670208.) - (110934582115.*PI**2.*eta/512004096.) - (1484623162301215.*eta**2./6604468835328.) + (128895671353745.*eta**3./217729741824.) \
	    - (1140350944327.*log(2.)/24000192.) + (1296725746149.*log(3.)/49778176.))*vByv0**(-37./3.) + ((-58408920984602265817973./7199481785765068800.) + (300051120571.*PI**2./970776576.) + (211649317.*GammaE/191520.) \
	    - (321242351927990989.*eta/133571090644992.) + (6450375061.*PI**2.*eta/26966016.) + (121407711933667.*eta**2./144557457408.) + (49171400252465.*eta**3./91738386432.) + (2117998887803.*log(2.)/44241120.) - (334711679031.*log(3.)/13108480.))*vByv0**(-38./3.) \
	    + ((1186114296954056489./17424955915960320.) + (505927225190405411.*eta/622319854141440.) - (2049675663640763.*eta**2./2469523230720.) - (40063118477671.*eta**3./20353213440.))*vByv0**(-44./3.) + (-2341612230425.*PI**2/3675082752.)*vByv0**(-47./3.) \
	    + ((4305919023475945625./31959919143419904.) - (635845376630254675.*eta/1141425683693568.) - (1251924162416035.*eta**2./1509822332928.) + (91862546967565.*eta**3./37330771968.))*vByv0**(-50./3.) + ((259620437372696563./159257838845952.) \
	    + (691917129965.*PI**2./2589262848.) - (558835855.*GammaE/2030112.) - (245999063921173.*eta/13702378991616.) - (20770936405.*PI**2.*eta/575391744.) + (255806950720535.*eta**2./326247118848.) - (9022269087085.*eta**3./8738762112.) - (12629690323.*log(2.)/188800416.) \
	    - (27159422553.*log(3.)/55940864.))*vByv0**(-56./3.))*e0**4.) + ((((20158674516353278980289./4685262871855104000.) - (7268851140841.*PI**2./11498029056.) - (248175681337.*GammaE/698664960.) - (691570196940108095.*eta/121695139528704.) + (174242180275.*PI**2.*eta/1277558784.) \
	    - (1597417452214177.*eta**2./1448751661056.) - (116996625810085.*eta**3./155223392256.) - (3265989073483.*log(2.)/2095994880.) + (172815325821.*log(3.)/207011840.))*vByv0**(-19./3.) + (-(345489155963130081415./616542346892279808.) - (2226563286362925308965.*eta/1078949107061489664.) \
	    + (6192425181604329515.*eta**2./4281544075640832.) + (129063292052563975.*eta**3./35287451172864.))*vByv0**(-25./3.) + (10562258458043923.*PI**2./7085660405760.)*vByv0**(-28./3.) + (-(1008593585234921446306165./4176146874611747782656.) + (54557236190663187381245.*eta/21306871809243611136.) \
	    + (3698038495111948565.*eta**2./2167976374566912.) - (30328195477605980725.*eta**3./4877946842775552.))*vByv0**(-31./3.) + ((1720927919854684009084595897./40516888294827538513920.) - (3986831179520597.*PI**2./8405059239936.) - (1423526912698421.*GammaE/255362042880.) \
	    + (131543151853096063653535.*eta/6190510052685643776.) - (775201866281389.*PI**2.*eta/466947735552.) - (174156473319672237061.*eta**2./96372409245106176.) + (18782995537481836405.*eta**3./5162807638130688.) - (233279096767651103.*log(2.)/766086128640.) + (773344207011339.*log(3.)/4450754560.) \
            - (164052734375.*log(5.)/20959232.))*vByv0**(-37./3.) - (194092844431833329313124279./3282963694308871372800.) + ((997069873657433.*PI**2./442674118656.) - ((1067488335456714056447.*eta)/60908417334116352.) + ((403437826755575441.*eta**2.)/65918200578048.) + ((163396563038941195.*eta**3.)/41832704212992.) \
	    + ((703310680391.*GammaE)/87333120.) + ((21434596327703.*eta*PI**2.)/12296503296.) + ((7038110304169369.*log(2.))/20173950720.) - ((370748969806671.*log(3.))/1992488960.))*vByv0**(-38./3.) + ((157835445312097689526421./206590277339625553920.) + ((68920046343599181068441.*eta)/7378224190700912640.) \
	    - ((1578991773134339939.*eta**2.)/250245020712960.) - ((350432910788522809.*eta**3.)/18562130657280.))*vByv0**(-44./3.) + (-(39680793155110375.*PI**2./6703350939648.))*vByv0**(-47./3.) + ((1670461618136994864075875./663104402387676168192.) - ((172890434755684698506065.*eta)/23682300085274148864.) \
	    - ((8402731601655384835.*eta**2.)/556079179235328.) + ((1860567315439539235.*eta**3.)/59579912060928.))*vByv0**(-50./3.) + (-(95765636723679036324982133./3502538538798360821760.) - ((55011254544918787424693.*eta)/2274375674544390144.) + ((463380478491174152645.*eta**2.)/27075900887433216.) \
	    - ((3879443939044136875.*eta**3.)/223153029292032.) + ((45970619802497.*GammaE)/14348831616.) + ((13756565834952955.*PI**2.)/4722815434752.) + ((20907767625235.*eta*PI**2.)/16398664704.) + ((13266735591208763.*log(2.))/43046494848.) - ((726469287588495.*log(3.))/4251505664.))*vByv0**(-56./3.) \
	    + ((480248792373794040814145387./9311678842039707893760.) + ((40273982144768174907437.*eta)/2015514900874395648.) - ((404897073357722903.*eta**2.)/68359615414272.) - ((3388909956719855.*eta**3.)/813804945408.) - ((5813865129161.*GammaE)/815109120.) - ((1918746213416491.*PI**2.)/1000328527872.)\
	    - ((7551348704987.*eta*PI**2.)/4471455744.) - ((22635681300089561.*log(2.))/66023838720.) + ((18242991444087.*log(3.))/103505920.) + ((910126953125.*log(5.))/146313216.))*vByv0**(-19.) + (-(2789890631138467908605./13735360805958844416.) - ((3754327126582174511645.*eta)/490548600212815872.) \
	    + ((872784891185699585.*eta**2.)/149740109955072.) + ((34162866264153035.*eta**3.)/1782620356608.))*vByv0**(-21.) + (398174549166095.*PI**2./80486203392.)*vByv0**(-22.) + (-(5758325856259741./2983951466496.) + ((1004822297691241255.*eta)/214844505587712.) + ((12453157351854917.*eta**2.)/852557561856.) \
	    - ((35025987744365.*eta**3.)/1171095552.))*vByv0**(-23.) + (-(1844247076182880525./167330816851968.) + ((2728185267249325.*eta)/633828851712.) - ((1304478350387875.*eta**2.)/80486203392.) + ((134782341955625.*eta**3.)/8623521792.) + ((249956266625.*GammaE)/139732992.) - ((10225600094125.*PI**2.)/3832676352.) \
	    + ((299691309125.*eta*PI**2.)/1277558784.) + ((182226181475.*log(2.))/419198976.) + ((130622307075.*log(3.))/41402368.))*vByv0**(-25.))*e0**6.)
      
    p6Lecc =  (((-734341./16800.)*log(v)*vByv0**(-19./3.) + (2603845./61404.)*(log(v) - log(vByv0))*vByv0**(-37./3.))*e0**2.) + (((-2440215143./15321600.)*vByv0**(-19./3.)*log(v) + (36290762107./56000448.)*(log(vByv0) - log(v))*vByv0**(-37./3.) \
	      + (211649317./191520.)*log(v)*vByv0**(-38./3.) + (558835855./2030112.)*(log(vByv0) - log(v))*vByv0**(-56./3.))*e0**4.) + ((-(248175681337./698664960.)*log(v)*vByv0**(-19./3.) + (1423526912698421./255362042880.)*(log(vByv0) - log(v))*vByv0**(-37./3.) \
              + (703310680391./87333120.)*log(v)*vByv0**(-38./3.) + (45970619802497./14348831616.)*(log(v) - log(vByv0))*vByv0**(-56./3.) - (5813865129161./815109120.)*log(v)*vByv0**(-19.) + (249956266625./139732992.)*(log(v) - log(vByv0))*vByv0**(-25.))*e0**6.)

'''
    #Circular + eccentric phasing    
    phase = 2*f*PI*tc - phic - PI/4. + (3./(128.*v**5*eta))*(p0 + p0ecc + v*(p1 + p1ecc) + v**2*(p2 + p2ecc) + v**3*(p3 + p3ecc) + v**4*(p4 + p4ecc) + v**5*(p5 + p5L + p5ecc + p5Lecc) + v**6*(p6 + p6L + p6ecc + p6Lecc) + v**7*p7)
'''  
    #phase due to tidal heating
    #----------------------------------------------------
    psi_so1 = (1/6.)*(-56*eta - 73*np.sqrt(1 - 4*eta) + 73)*chi1z
    psi_so2 = (1/6.)*(-56*eta + 73*np.sqrt(1 - 4*eta) + 73)*chi2z
    psi_so = psi_so1 + psi_so2
    #con = (3./(128.*v**5*eta))
    THterm1 = -(10/9.)*(v**5)*Heff5*(3*np.log(v))
    THterm2 = -(5/168.)*(v**7)*Heff5*(952*eta + 995)
    THterm3 = (5/9.)*(v**8)*(3*np.log(v))*(-4*Heff8 + Heff5*psi_so)
    
    #heated_phase = con*(term1 + term2 + term3)
    
    #Circular + eccentric phasing with Tidal Heating terms at 2.5, 3.5 & 4PN orders    
    phase1 = 2*f*PI*tc - 1*phic - PI/4. + (3./(128.*v1**5*eta))*(p0 + p0ecc + v1*(p1 + p1ecc) + v1**2*(p2 + p2ecc) + v1**3*(p3 + p3ecc) + v1**4*(p4 + p4ecc) + v1**5*(p5 + p5L + p5ecc + p5Lecc + THterm1) + v1**6*(p6 + p6L + p6ecc + p6Lecc) + v1**7*(p7 + THterm2) + v1**8*THterm3)
    phase2 = 2*f*PI*tc - 2*phic - PI/4. + (3./(128.*v2**5*eta))*(p0 + p0ecc + v2*(p1 + p1ecc) + v2**2*(p2 + p2ecc) + v2**3*(p3 + p3ecc) + v2**4*(p4 + p4ecc) + v2**5*(p5 + p5L + p5ecc + p5Lecc + THterm1) + v2**6*(p6 + p6L + p6ecc + p6Lecc) + v2**7*(p7 + THterm2) + v2**8*THterm3)
    phase3 = 2*f*PI*tc - 3*phic - PI/4. + (3./(128.*v3**5*eta))*(p0 + p0ecc + v3*(p1 + p1ecc) + v3**2*(p2 + p2ecc) + v3**3*(p3 + p3ecc) + v3**4*(p4 + p4ecc) + v3**5*(p5 + p5L + p5ecc + p5Lecc + THterm1) + v3**6*(p6 + p6L + p6ecc + p6Lecc) + v3**7*(p7 + THterm2) + v3**8*THterm3)
    phase4 = 2*f*PI*tc - 4*phic - PI/4. + (3./(128.*v4**5*eta))*(p0 + p0ecc + v4*(p1 + p1ecc) + v4**2*(p2 + p2ecc) + v4**3*(p3 + p3ecc) + v4**4*(p4 + p4ecc) + v4**5*(p5 + p5L + p5ecc + p5Lecc + THterm1) + v4**6*(p6 + p6L + p6ecc + p6Lecc) + v4**7*(p7 + THterm2) + v4**8*THterm3)
    phase5 = 2*f*PI*tc - 5*phic - PI/4. + (3./(128.*v5**5*eta))*(p0 + p0ecc + v5*(p1 + p1ecc) + v5**2*(p2 + p2ecc) + v5**3*(p3 + p3ecc) + v5**4*(p4 + p4ecc) + v5**5*(p5 + p5L + p5ecc + p5Lecc + THterm1) + v5**6*(p6 + p6L + p6ecc + p6Lecc) + v5**7*(p7 + THterm2) + v5**8*THterm3)
    phase6 = 2*f*PI*tc - 6*phic - PI/4. + (3./(128.*v6**5*eta))*(p0 + p0ecc + v6*(p1 + p1ecc) + v6**2*(p2 + p2ecc) + v6**3*(p3 + p3ecc) + v6**4*(p4 + p4ecc) + v6**5*(p5 + p5L + p5ecc + p5Lecc + THterm1) + v6**6*(p6 + p6L + p6ecc + p6Lecc) + v6**7*(p7 + THterm2) + v6**8*THterm3)
    phase7 = 2*f*PI*tc - 7*phic - PI/4. + (3./(128.*v7**5*eta))*(p0 + p0ecc + v7*(p1 + p1ecc) + v7**2*(p2 + p2ecc) + v7**3*(p3 + p3ecc) + v7**4*(p4 + p4ecc) + v7**5*(p5 + p5L + p5ecc + p5Lecc + THterm1) + v7**6*(p6 + p6L + p6ecc + p6Lecc) + v7**7*(p7 + THterm2) + v7**8*THterm3)
    phase8 = 2*f*PI*tc - 8*phic - PI/4. + (3./(128.*v8**5*eta))*(p0 + p0ecc + v8*(p1 + p1ecc) + v8**2*(p2 + p2ecc) + v8**3*(p3 + p3ecc) + v8**4*(p4 + p4ecc) + v8**5*(p5 + p5L + p5ecc + p5Lecc + THterm1) + v8**6*(p6 + p6L + p6ecc + p6Lecc) + v8**7*(p7 + THterm2) + v8**8*THterm3)

    #phase += heated_phase
    
    hp = A * f**(-7./6.) * (hp1 * (1/2)**(2./3.) * (cos(phase1) + 1j*sin(phase1)) * ustep(1*flso - 2*f) + hp2 * (2/2)**(2./3.) * (cos(phase2) + 1j*sin(phase2)) * ustep(2*flso - 2*f) + hp3 * (3/2)**(2./3.) * (cos(phase3) + 1j*sin(phase3)) * ustep(3*flso - 2*f) + hp4 * (4/2)**(2./3.) * (cos(phase4) + 1j*sin(phase4)) * ustep(4*flso - 2*f) + hp5 * (5/2)**(2./3.) * (cos(phase5) + 1j*sin(phase5)) * ustep(5*flso - 2*f) + hp6 * (6/2)**(2./3.) * (cos(phase6) + 1j*sin(phase6)) * ustep(6*flso - 2*f) + hp7 * (7/2)**(2./3.) * (cos(phase7) + 1j*sin(phase7)) * ustep(7*flso - 2*f) + hp8 * (8/2)**(2./3.) * (cos(phase8) + 1j*sin(phase8)) * ustep(8*flso - 2*f))	 #Plus GW polarization for eccentric binaries 
    hc = A * f**(-7./6.) * (hc1 * (1/2)**(2./3.) * (cos(phase1) + 1j*sin(phase1)) * ustep(1*flso - 2*f) + hc2 * (2/2)**(2./3.) * (cos(phase2) + 1j*sin(phase2)) * ustep(2*flso - 2*f) + hc3 * (3/2)**(2./3.) * (cos(phase3) + 1j*sin(phase3)) * ustep(3*flso - 2*f) + hc4 * (4/2)**(2./3.) * (cos(phase4) + 1j*sin(phase4)) * ustep(4*flso - 2*f) + hc5 * (5/2)**(2./3.) * (cos(phase5) + 1j*sin(phase5)) * ustep(5*flso - 2*f) + hc6 * (6/2)**(2./3.) * (cos(phase6) + 1j*sin(phase6)) * ustep(6*flso - 2*f) + hc7 * (7/2)**(2./3.) * (cos(phase7) + 1j*sin(phase7)) * ustep(7*flso - 2*f) + hc8 * (8/2)**(2./3.) * (cos(phase8) + 1j*sin(phase8)) * ustep(8*flso - 2*f))     #Cross GW polarization for eccentric binaries
	
  #  hp = 0.5*(1+(cos(iota))**2)*A*f**(-7./6.)*(cos(phase) - 1j*sin(phase))    Plus GW polarization for circular binaries 
  #  hc = -1j*cos(iota)*A*f**(-7./6.)*(cos(phase) - 1j*sin(phase))             Cross GW polarization for circular binaries

    return hp, hc
