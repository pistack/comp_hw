/*!
 * @file node_weight_table.hpp
 * @ingroup libpath
 * @brief table for node and weights
 * @author pistack (Junho Lee)
 * @date 2021. 11. 10.
 */

#ifndef NODE_WEIGHT_TABLE_H
#define NODE_WEIGHT_TABLE_H
#include <vector>

namespace libpath{
/// @brief table for gauss kronrod node and weights
/// @param T precision should be one of float, double, long double
/// @param N order of gauss-kronrod quadrature
/// currently only supports N=15, 21, 31, 41, 51, 61
/// @note coefficients are obtained from 
/// [gau-kronrod-nodes-weights]
/// (https://www.advanpix.com/2011/11/07/gauss-kronrod-quadrature-nodes-weights/)
/// @ingroup libpath
template<typename T, int N>
class gau_kron_table
{
    public:
    int order; ///< order of gauss-kronrod quadrature;
    std::vector<T> nodes; ///< node of gauss-kronrod quadrature;
    std::vector<T> weight_gauss; ///< weight of gauss quadrature;
    std::vector<T> weight_kronrod; ///< weight of kronrod quadrature;
};

#ifndef DOXYGEN_SKIP
template<typename T>
class gau_kron_table<T, 15>
{
    public:
    const int order = 15;
    const std::vector<T> nodes
    {
        9.914553711208126392068546975263285e-01l,
        9.491079123427585245261896840478513e-01l,
        8.648644233597690727897127886409262e-01l,
        7.415311855993944398638647732807884e-01l,
        5.860872354676911302941448382587296e-01l,
        4.058451513773971669066064120769615e-01l,
        2.077849550078984676006894037732449e-01l
    };
    const std::vector<T> weight_gauss
    {
        1.294849661688696932706114326790820e-01l,
        2.797053914892766679014677714237796e-01l,
        3.818300505051189449503697754889751e-01l,
        4.179591836734693877551020408163265e-01l
    };
    const std::vector<T> weight_kronrod
    {
        2.293532201052922496373200805896959e-02l,
        6.309209262997855329070066318920429e-02l,
        1.047900103222501838398763225415180e-01l,
        1.406532597155259187451895905102379e-01l,
        1.690047266392679028265834265985503e-01l,
        1.903505780647854099132564024210137e-01l,
        2.044329400752988924141619992346491e-01l,
        2.094821410847278280129991748917143e-01l
    };
};

template<typename T>
class gau_kron_table<T, 21>
{
    public:
    const int order = 21;
    const std::vector<T> nodes
    {
        9.956571630258080807355272806890028e-01l,
        9.739065285171717200779640120844521e-01l,
        9.301574913557082260012071800595083e-01l,
        8.650633666889845107320966884234930e-01l,
        7.808177265864168970637175783450424e-01l,
        6.794095682990244062343273651148736e-01l,
        5.627571346686046833390000992726941e-01l,
        4.333953941292471907992659431657842e-01l,
        2.943928627014601981311266031038656e-01l,
        1.488743389816312108848260011297200e-01l
    };
    const std::vector<T> weight_gauss
    {
        6.667134430868813759356880989333179e-02l,
        1.494513491505805931457763396576973e-01l,
        2.190863625159820439955349342281632e-01l,
        2.692667193099963550912269215694694e-01l,
        2.955242247147528701738929946513383e-01l
    };
    const std::vector<T> weight_kronrod
    {
        1.169463886737187427806439606219205e-02l,
        3.255816230796472747881897245938976e-02l,
        5.475589657435199603138130024458018e-02l,
        7.503967481091995276704314091619001e-02l,
        9.312545458369760553506546508336634e-02l,
        1.093871588022976418992105903258050e-01l,
        1.234919762620658510779581098310742e-01l,
        1.347092173114733259280540017717068e-01l,
        1.427759385770600807970942731387171e-01l,
        1.477391049013384913748415159720680e-01l,
        1.494455540029169056649364683898212e-01l
    };
};

template<typename T>
class gau_kron_table<T, 31>
{
    public:
    const int order = 31;
    const std::vector<T> nodes
    {
        9.980022986933970602851728401522712e-01l, 
        9.879925180204854284895657185866126e-01l,
        9.677390756791391342573479787843372e-01l, 
        9.372733924007059043077589477102095e-01l,
        8.972645323440819008825096564544959e-01l, 
        8.482065834104272162006483207742169e-01l,
        7.904185014424659329676492948179473e-01l, 
        7.244177313601700474161860546139380e-01l,
        6.509967412974169705337358953132747e-01l, 
        5.709721726085388475372267372539106e-01l,
        4.850818636402396806936557402323506e-01l, 
        3.941513470775633698972073709810455e-01l,
        2.991800071531688121667800242663890e-01l, 
        2.011940939974345223006283033945962e-01l,
        1.011420669187174990270742314473923e-01l	
    };
    const std::vector<T> weight_gauss
    {
        3.075324199611726835462839357720442e-02l, 
        7.036604748810812470926741645066734e-02l,
        1.071592204671719350118695466858693e-01l, 
        1.395706779261543144478047945110283e-01l,
        1.662692058169939335532008604812088e-01l, 
        1.861610000155622110268005618664228e-01l,
        1.984314853271115764561183264438393e-01l, 
        2.025782419255612728806201999675193e-01l
    };
    const std::vector<T> weight_kronrod
    {
        5.377479872923348987792051430127650e-03l, 
        1.500794732931612253837476307580727e-02l,
        2.546084732671532018687400101965336e-02l, 
        3.534636079137584622203794847836005e-02l,
        4.458975132476487660822729937327969e-02l, 
        5.348152469092808726534314723943030e-02l,
        6.200956780067064028513923096080293e-02l, 
        6.985412131872825870952007709914748e-02l,
        7.684968075772037889443277748265901e-02l, 
        8.308050282313302103828924728610379e-02l,
        8.856444305621177064727544369377430e-02l, 
        9.312659817082532122548687274734572e-02l,
        9.664272698362367850517990762758934e-02l, 
        9.917359872179195933239317348460313e-02l,
        1.007698455238755950449466626175697e-01l, 
        1.013300070147915490173747927674925e-01l
    };
};

template<typename T>
class gau_kron_table<T, 41>
{
    public:
    const int order=41;
    const std::vector<T> nodes
    {
        9.988590315882776638383155765458630e-01l,
        9.931285991850949247861223884713203e-01l,
        9.815078774502502591933429947202169e-01l,
        9.639719272779137912676661311972772e-01l,
        9.408226338317547535199827222124434e-01l,
        9.122344282513259058677524412032981e-01l,
        8.782768112522819760774429951130785e-01l,
        8.391169718222188233945290617015207e-01l,
        7.950414288375511983506388332727879e-01l,
        7.463319064601507926143050703556416e-01l,
        6.932376563347513848054907118459315e-01l,
        6.360536807265150254528366962262859e-01l,
        5.751404468197103153429460365864251e-01l,
        5.108670019508270980043640509552510e-01l,
        4.435931752387251031999922134926401e-01l,
        3.737060887154195606725481770249272e-01l,
        3.016278681149130043205553568585923e-01l,
        2.277858511416450780804961953685746e-01l,
        1.526054652409226755052202410226775e-01l,
        7.652652113349733375464040939883821e-02l	
    };
    const std::vector<T> weight_gauss
    {
        1.761400713915211831186196235185282e-02l,
        4.060142980038694133103995227493211e-02l,
        6.267204833410906356950653518704161e-02l,
        8.327674157670474872475814322204621e-02l,
        1.019301198172404350367501354803499e-01l,
        1.181945319615184173123773777113823e-01l,
        1.316886384491766268984944997481631e-01l,
        1.420961093183820513292983250671649e-01l,
        1.491729864726037467878287370019694e-01l,
        1.527533871307258506980843319550976e-01l
    };
    const std::vector<T> weight_kronrod
    {
        3.073583718520531501218293246030987e-03l,
        8.600269855642942198661787950102347e-03l,
        1.462616925697125298378796030886836e-02l,
        2.038837346126652359801023143275471e-02l,
        2.588213360495115883450506709615314e-02l,
        3.128730677703279895854311932380074e-02l,
        3.660016975820079803055724070721101e-02l,
        4.166887332797368626378830593689474e-02l,
        4.643482186749767472023188092610752e-02l,
        5.094457392372869193270767005034495e-02l,
        5.519510534828599474483237241977733e-02l,
        5.911140088063957237496722064859422e-02l,
        6.265323755478116802587012217425498e-02l,
        6.583459713361842211156355696939794e-02l,
        6.864867292852161934562341188536780e-02l,
        7.105442355344406830579036172321017e-02l,
        7.303069033278666749518941765891311e-02l,
        7.458287540049918898658141836248753e-02l,
        7.570449768455667465954277537661656e-02l,
        7.637786767208073670550283503806100e-02l,
        7.660071191799965644504990153010174e-02l
    };
};

template<typename T>
class gau_kron_table<T, 51>
{
    public:
    const int order=51;
    const std::vector<T> nodes
    {
        9.992621049926098341934574865403406e-01l,
        9.955569697904980979087849468939016e-01l,
        9.880357945340772476373310145774062e-01l,
        9.766639214595175114983153864795941e-01l,
        9.616149864258425124181300336601672e-01l,
        9.429745712289743394140111696584705e-01l,
        9.207471152817015617463460845463306e-01l,
        8.949919978782753688510420067828050e-01l,
        8.658470652932755954489969695883401e-01l,
        8.334426287608340014210211086935696e-01l,
        7.978737979985000594104109049943066e-01l,
        7.592592630373576305772828652043610e-01l,
        7.177664068130843881866540797732978e-01l,
        6.735663684734683644851206332476222e-01l,
        6.268100990103174127881226816245179e-01l,
        5.776629302412229677236898416126541e-01l,
        5.263252843347191825996237781580102e-01l,
        4.730027314457149605221821150091920e-01l,
        4.178853821930377488518143945945725e-01l,
        3.611723058093878377358217301276407e-01l,
        3.030895389311078301674789099803393e-01l,
        2.438668837209884320451903627974516e-01l,
        1.837189394210488920159698887595284e-01l,
        1.228646926107103963873598188080368e-01l,
        6.154448300568507888654639236679663e-02l	
    };
    const std::vector<T> weight_gauss
    {
        1.139379850102628794790296411323477e-02l,
        2.635498661503213726190181529529914e-02l,
        4.093915670130631265562348771164595e-02l,
        5.490469597583519192593689154047332e-02l,
        6.803833381235691720718718565670797e-02l,
        8.014070033500101801323495966911130e-02l,
        9.102826198296364981149722070289165e-02l,
        1.005359490670506442022068903926858e-01l,
        1.085196244742636531160939570501166e-01l,
        1.148582591457116483393255458695558e-01l,
        1.194557635357847722281781265129010e-01l,
        1.222424429903100416889595189458515e-01l,
        1.231760537267154512039028730790501e-01l
    };
    const std::vector<T> weight_kronrod
    {
        1.987383892330315926507851882843410e-03l,
        5.561932135356713758040236901065522e-03l,
        9.473973386174151607207710523655324e-03l,
        1.323622919557167481365640584697624e-02l,
        1.684781770912829823151666753633632e-02l,
        2.043537114588283545656829223593897e-02l,
        2.400994560695321622009248916488108e-02l,
        2.747531758785173780294845551781108e-02l,
        3.079230016738748889110902021522859e-02l,
        3.400213027432933783674879522955120e-02l,
        3.711627148341554356033062536761988e-02l,
        4.008382550403238207483928446707565e-02l,
        4.287284502017004947689579243949516e-02l,
        4.550291304992178890987058475266039e-02l,
        4.798253713883671390639225575691475e-02l,
        5.027767908071567196332525943344008e-02l,
        5.236288580640747586436671213787271e-02l,
        5.425112988854549014454337045987561e-02l,
        5.595081122041231730824068638274735e-02l,
        5.743711636156783285358269393950647e-02l,
        5.868968002239420796197417585678776e-02l,
        5.972034032417405997909929193256185e-02l,
        6.053945537604586294536026751756543e-02l,
        6.112850971705304830585903041629271e-02l,
        6.147118987142531666154413196526418e-02l,
        6.158081806783293507875982424006455e-02l
    };
};

template<typename T>
class gau_kron_table<T, 61>
{
    public:
    const int order = 61;
    const std::vector<T> nodes
    {
        9.994844100504906375713258957058108e-01l,
        9.968934840746495402716300509186953e-01l,
        9.916309968704045948586283661094857e-01l,
        9.836681232797472099700325816056628e-01l,
        9.731163225011262683746938684237069e-01l,
        9.600218649683075122168710255817977e-01l,
        9.443744447485599794158313240374391e-01l,
        9.262000474292743258793242770804740e-01l,
        9.055733076999077985465225589259583e-01l,
        8.825605357920526815431164625302256e-01l,
        8.572052335460610989586585106589439e-01l,
        8.295657623827683974428981197325019e-01l,
        7.997278358218390830136689423226832e-01l,
        7.677774321048261949179773409745031e-01l,
        7.337900624532268047261711313695276e-01l,
        6.978504947933157969322923880266401e-01l,
        6.600610641266269613700536681492708e-01l,
        6.205261829892428611404775564311893e-01l,
        5.793452358263616917560249321725405e-01l,
        5.366241481420198992641697933110728e-01l,
        4.924804678617785749936930612077088e-01l,
        4.470337695380891767806099003228540e-01l,
        4.004012548303943925354762115426606e-01l,
        3.527047255308781134710372070893739e-01l,
        3.040732022736250773726771071992566e-01l,
        2.546369261678898464398051298178051e-01l,
        2.045251166823098914389576710020247e-01l,
        1.538699136085835469637946727432559e-01l,
        1.028069379667370301470967513180006e-01l,
        5.147184255531769583302521316672257e-02l
    };
    const std::vector<T> weight_gauss
    {
        7.968192496166605615465883474673622e-03,
        1.846646831109095914230213191204727e-02,
        2.878470788332336934971917961129204e-02,
        3.879919256962704959680193644634769e-02,
        4.840267283059405290293814042280752e-02,
        5.749315621761906648172168940205613e-02,
        6.597422988218049512812851511596236e-02,
        7.375597473770520626824385002219073e-02,
        8.075589522942021535469493846052973e-02,
        8.689978720108297980238753071512570e-02,
        9.212252223778612871763270708761877e-02,
        9.636873717464425963946862635180987e-02,
        9.959342058679526706278028210356948e-02,
        1.017623897484055045964289521685540e-01,
        1.028526528935588403412856367054150e-01
    };
    const std::vector<T> weight_kronrod
    {
        1.389013698677007624551591226759700e-03l,
        3.890461127099884051267201844515503e-03l,
        6.630703915931292173319826369750168e-03l,
        9.273279659517763428441146892024360e-03l,
        1.182301525349634174223289885325059e-02l,
        1.436972950704580481245143244358001e-02l,
        1.692088918905327262757228942032209e-02l,
        1.941414119394238117340895105012846e-02l,
        2.182803582160919229716748573833899e-02l,
        2.419116207808060136568637072523203e-02l,
        2.650995488233310161060170933507541e-02l,
        2.875404876504129284397878535433421e-02l,
        3.090725756238776247288425294309227e-02l,
        3.298144705748372603181419101685393e-02l,
        3.497933802806002413749967073146788e-02l,
        3.688236465182122922391106561713597e-02l,
        3.867894562472759295034865153228105e-02l,
        4.037453895153595911199527975246811e-02l,
        4.196981021516424614714754128596976e-02l,
        4.345253970135606931683172811707326e-02l,
        4.481480013316266319235555161672324e-02l,
        4.605923827100698811627173555937358e-02l,
        4.718554656929915394526147818109949e-02l,
        4.818586175708712914077949229830459e-02l,
        4.905543455502977888752816536723817e-02l,
        4.979568342707420635781156937994233e-02l,
        5.040592140278234684089308565358503e-02l,
        5.088179589874960649229747304980469e-02l,
        5.122154784925877217065628260494421e-02l,
        5.142612853745902593386287921578126e-02l,
        5.149472942945156755834043364709931e-02l
    };
};
#endif /* DOXYGEN_SKIP */
}

#endif