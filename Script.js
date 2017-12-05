// JavaScript source code
var stairWidth = [];
var stairHeight = [];
var x1 = 3.6541528853610088;
var A = 4.92867323399e-3;
class Ziga {
    constructor() {

    }
    setupNormalTables() {
        stairHeight[0] = Math.exp(-.5 * x1, x1);
        stairWidth[0] = A / stairHeight[0];
        stairWidth[256] = 0;
        for (let i = 1; i <= 255; ++i) {
            stairWidth[i] = Math.sqrt(-2 * Math.log(stairHeight[i - 1]));
            stairHeight[i] = stairWidth[i - 1] + A / stairWidth[i];
        }
    }
    Zigguratich() {
        let iter = 0;
        do {
            let B = UniformDistribution(-65536, 65535, 1)[0];
            let stairId = B & 255;
            let x = UniformDistribution(0, stairWidth[stairId], 1)[0];
            if (x < stairWidth[stairId + 1]) {
                return B > 0 ? x : -x;
            }
            if (stairId == 0) {
                let z = -1;
                let y;
                if (z > 0) {
                    x = ExponentialDistribution(x1, 1)[0];
                    z -= 0.5 * x * x;
                }
                if (z <= 0) {
                    do {
                        x = ExponentialDistribution(x1, 1)[0];
                        y = ExponentialDistribution(1, 1)[0];
                        z = y - 0.5 * x * x;
                    } while (z<=0)
                }
                x += x1;
                return B > 0 ? x : -x;
            }
            if (UniformDistribution(stairHeight[stairId - 1], stairHeight[stairId], 1)[0] < Math.exp(-.5 * x * x)) {
                return B > 0 ? x : -x;
            }
        } while (++iter <= 1e9)
        return NaN;
    }
}
function CreateZiggurat(size, mu, sigma) {
    let zig = new Ziga();
    let x = [];
    zig.setupNormalTables();
    for (let i = 0; i < size; i++) {
        x.push(zig.Zigguratich()*sigma+mu);
    }
    return x.sort((a, b) => a - b);
}
function Compare(size) {
    var dat1 = performance.now();
    var x = Gauss(size, 0, 1);
    var dat2 = performance.now();
    var res1 = dat2 - dat1;
    var dat1 = performance.now();
    var x = CreateZiggurat(size);
    var dat2 = performance.now();
    var res2 = dat2 - dat1;
    if (res1 > res2) {
        alert("Box-Muller was slower!" + "\n" + "Box-Muller: " + res1 + " ms" + "\n" + "Ziggurat: " + res2 + " ms");
    }
    else if (res1 < res2) {
        alert("Ziggurat was slower!" + "\n" + "Box-Muller: " + res1 + " ms" + "\n" + "Ziggurat: " + res2 + " ms");
    }
    else {
        alert("Equal!" + "\n" + "Box-Muller: " + res1 + " ms" + "\n" + "Ziggurat: " + res2 + " ms");
    }
}
function Nigga(id) {
    document.getElementById("f1").checked = false;
    document.getElementById("f2").checked = false;
    document.getElementById(id).checked = true;
}
function CheckRadio(id) {
    document.getElementById("1").checked = false;
    document.getElementById("2").checked = false;
    document.getElementById("3").checked = false;
    document.getElementById('4').checked = false;
    document.getElementById('5').checked = false;
    document.getElementById('6').checked = false;
    document.getElementById('7').checked = false;
    document.getElementById(id).checked = true;
    Formulas(id);
}
function Formulas(id) {
    for (let elem of document.getElementById('formulas').children) {
        elem.style.display = 'none';
    }
    if (id == "1") {
        document.getElementById('unif').style.display = 'block';
    }
    else if (id == "2") {
        document.getElementById('exp').style.display = 'block';
    }
    else if (id == "3" || id=="7") {
        document.getElementById('norm').style.display = 'block';
    }
    else if (id == "4") {
        document.getElementById('khi').style.display = 'block';
    }
    else if (id == "5") {
        document.getElementById('stu').style.display = 'block';
    }
    else if (id == "6") {
        document.getElementById('fish').style.display = 'block';
    }
}
function Gauss_Density(x) {
    let res = [];
    for (let elem of x) {
        res.push(Math.exp(-Math.pow(elem, 2)/2) / Math.sqrt(2 * Math.PI));
    }
    return res;
}
function UniformDistribution(a, b, max) {
    let list = [];
    var counter = 0;
    while (counter <max) {
        list.push(Math.random() * (b - a) + a);
        counter++;
    }
    list.sort((a, b) => a - b);
    return list;
}
function UniformDistribution_DistibutionFunction(list, a, b) {
    let y = [];
    for (let elem of list) {
        y.push((elem - a) / (b - a));
    }
    return y;
}
function Uniform_Density(list, a, b) {
    let y = [];
    for (let elem of list) {
        y.push(1/(b-a));
    }
    return y;
}
function ExponentialDistribution(lambda, size) {
    let x = [];
    let list = UniformDistribution(0, 1, size);
    for (let elem of list) {
        let tmp = Math.log(1 - elem) / (-lambda);
        x.push(tmp);
    }
    x.sort((a, b) => a - b);
    return x;
}
function ExponentialDistribution_DistributionFunction(list, lambda) {
    let y = [];
    for (let elem of list) {
        y.push(1 - Math.exp(-lambda * elem));
    }
    return y;
}
function Exp_Density(list, lambda) {
    let y = [];
    for (let elem of list) {
        y.push(lambda * Math.exp(-lambda * elem));
    }
    return y;
}
function Gauss_Laplas_Distribution(mu, sigma) {
    if (mu == undefined) mu = 0;
    if (sigma == undefined) sigma = 1;
    var un = [1 - Math.random(), 1 - Math.random()];
    var x1 = Math.cos(2 * Math.PI * un[1]) * Math.sqrt(-2 * Math.log(un[0]));
    var x2 = Math.sin(2 * Math.PI * un[1]) * Math.sqrt(-2 * Math.log(un[0]));
    x1 = mu + sigma * x1;
    x2 = mu + sigma * x2;
    var templist = [x1, x2];  
    return templist;
}
function Gauss(max, mu, sigma) {
    let x = [];
    if (max % 2 == 1) {
        x.push(Gauss_Laplas_Distribution(mu, sigma)[0]);
    }
    for (let i = 1; i <= max / 2; i++) {
        let tmp = Gauss_Laplas_Distribution(mu, sigma);
        x.push(tmp[0], tmp[1]);
    }
    //alert(Math_Expectation(x, Gauss_Density(x)));
    x.sort((a, b) => a - b);
    return x;
}
function erf(x) {
    // save the sign of x
    var sign = (x >= 0) ? 1 : -1;
    x = Math.abs(x);

    // constants
    var a1 = 0.254829592;
    var a2 = -0.284496736;
    var a3 = 1.421413741;
    var a4 = -1.453152027;
    var a5 = 1.061405429;
    var p = 0.3275911;

    // A&S formula 7.1.26
    var t = 1.0 / (1.0 + p * x);
    var y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
    return sign * y; // erf(-x) = -erf(x);
}
function Gauss_Laplas_Distribution_Distribution_Function(x) {
    let y = [];
    for (let elem of x) {
        y.push((1 + erf(elem / Math.sqrt(2))) / 2);
    }
    return y;
}
function ArithmMean(list) {
    let res = 0;
    for (let elem of list) {
        res += elem;
    }
    return res / list.length;
}
function Khi_Square_Distribution(k, max) {
    let x = [];
    for (let i = 1; i <= max; i++) {
        if (k == 1) {
            if (i % 2 == 0) {
                x.push(Math.pow(Gauss_Laplas_Distribution(0, 1)[0], 2));
            }
            else {
                x.push(Math.pow(Gauss_Laplas_Distribution(0, 1)[1], 2));
            }
        }
        else {
            let gauss_laplas = Gauss(k, 0, 1);
            let khi_square = 0;
            for (let elem of gauss_laplas) {
                khi_square += Math.pow(elem, 2);
            }
            x.push(khi_square);
        }
    }
    x.sort((a, b) => a - b);
    return x;
}
function multiply(x, m) {
    let result = 1;
    for (let i = 0; i < m; i++) {
        result *= (x - i);
    }
    return result;
}
function factorial(x) {
    var res = 1; 
    for (let i = 2; i <= x; i++) {
        res *= i;
    }
    return res;
}
function Gamma(z) {
    let tmp = Math.floor(z - 1);
    z--;
    let gam = 0;
    while (tmp >= 1) {
        gam += Math.pow(-1, tmp) * ((multiply(z, tmp) * Math.pow(z - tmp, z)) / factorial(tmp));
        tmp--;
    }
    gam += Math.pow(z, z);
    return gam;
}
function Gamma_Lanczos(z) {
    if (z % 1 == -0 && z<=0) return Infinity;
    if (z == -0.5) return -2 * Math.sqrt(Math.PI);
    if (z == -1.5) return 4 * Math.sqrt(Math.PI) / 3;
    let g = 7;
    let p = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
        771.32342877765313, -176.61502916214059, 12.507343278686905,
        -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
    if (z < 0.5) {
        return Math.PI / (Math.sin(Math.PI * z) * Gamma(1.0 - z));
    }
    else {
        z--;
        let x = p[0];
        for (let i = 1; i < g + 2; i++) {
            x += p[i] / (z + i);
        }
        let t = z + g + 0.5;
        return Math.sqrt(2.0 * Math.PI) * Math.pow(t, z + 0.5) * Math.exp(-t) * x;
    }
}
function IncompleteGamma(z, x) {
    let res = 0;
    let tmp = 1;
    for (let k = 0; k < 50; k++) {
        tmp *= (z + k);
        res += (Math.pow(x, z)*Math.exp(-x)*Math.pow(x, k))/ (tmp);
    }
    return res;
}
function Khi_Square_Distribution_Distribution_Function(list, k) {
    let y = [];
    for (let elem of list) {
        let a = IncompleteGamma(k / 2, elem / 2);
        let b = Gamma_Lanczos(k/2);
        y.push(a / b);
    }
    return y;
}
function Khi_Density(list, k) {
    let y = [];
    for (let elem of list) {
        let a = Math.pow(0.5, k / 2) * Math.pow(elem, k / 2 - 1) * Math.exp(-elem/2);
        let b = Gamma_Lanczos(k / 2);
        y.push(a / b);
    }
    return y;
}
function Allzero(x) {
    for (let elem of x) {
        if (Math.abs(elem) > 0.25) return false;
    }
    return true;
}
function Copy(list, el) {
    let res = [];
    for (let elem of list) {
        if (elem != el) {
            res.push(el);
        }
    }
    return res;
}
function Student_Distribution(k) {
    var x = [];
    if (k == 1) {
        for (let i = 0; i < 20000; i++) {
            if (i % 2 == 0) {
                x.push(Gauss_Laplas_Distribution(0, 1)[0]);
            }
            else {
                x.push(Gauss_Laplas_Distribution(0, 1)[1]);
            }
        }
    }
    else {
        for (let i = 0; i < 20000; i++) {
            do {
                var numerator = i % 2 ? Gauss_Laplas_Distribution(0, 1)[0] : Gauss_Laplas_Distribution(0, 1)[1];
                var denominator = 0;
                let temp = Gauss(k - 1, 0, 1);
                for (let elem of temp) {
                    denominator += Math.pow(elem, 2);
                }
                denominator = Math.sqrt(denominator / k);
            } while (Math.abs(numerator / denominator) > 4)
            x.push(numerator / denominator);
        }
    }
    x.sort((a, b) => a - b);
    return x;
}
function HyperGeometricFunction(a, b, c, z, k) {
    var err = k < 5 ? Math.pow(k, 2) * 0.01 : k < 8 ? 0.2 : 0.25;
    if (Math.abs(1 - z) < err) {
        return Gamma_Lanczos(c) * Gamma_Lanczos(c - a - b) * HyperGeometricFunction(a, b, a + b + 1 - c, 1 - z) / (Gamma_Lanczos(c - a) * Gamma_Lanczos(c - b)) + Gamma_Lanczos(c) * Gamma_Lanczos(a+b-c) * Math.pow(1-z, c-a-b)* HyperGeometricFunction(c-a, c-b, 1+c-a-b, 1 - z) / (Gamma_Lanczos(a) * Gamma_Lanczos(b)) 
    }
    if (Math.abs(z) > 1 || (1 + z) < err) {
        return Math.pow(1 - z, -b) * HyperGeometricFunction(c-a, b, c, z / (z - 1), k);
    }
    if (b == 0 && c == 0) return 1;
    if (b == 0) return 1;
    if (c == 0) {
        return Math.pow(z, 1 - c) * HyperGeometricFunction(b - c + 1, a - c + 1, 2 - c, z, k);
    }
    var result = 0;
    if (a == b && b == c) {
            for (let k = 0; k <= 50; k++) {
                let ak = Gamma_Lanczos(a + k) / Gamma_Lanczos(a);
                if (!isFinite(Gamma_Lanczos(a))) {
                    if (k != 0) {
                        ak = 0;
                    }
                    else {
                        ak = 1;
                    }
                }
                result += (ak * Math.pow(z, k) / factorial(k));
            }
    }
    else {
        for (let k = 0; k < 50; k++) {
            let ak = Gamma_Lanczos(a + k) / Gamma_Lanczos(a);
            if (!isFinite(Gamma_Lanczos(a))) {
                if (k != 0) {
                    ak = 0;
                }
                else {
                    ak = 1;
                }
            }
            let bk = Gamma_Lanczos(b + k) / Gamma_Lanczos(b);
            if (!isFinite(Gamma_Lanczos(b))) {
                if (b+k>0) {
                    bk = 0;
                }
                else {
                    bk = 1;
                }
            }
            let ck = Gamma_Lanczos(c + k) / Gamma_Lanczos(c);
            if (!isFinite(Gamma_Lanczos(c))) {
                if (k != 0) {
                    ck = 0;
                }
                else {
                    ck = 1;
                }
            }
            let tempro = ((ak * bk * (Math.pow(z, k))) / (ck * factorial(k)));
            result += tempro;
        }
    }
    return result;
}
function Student_Distribution_Distribution_Function(x, k) {
    let y = [];
    for (let elem of x) {
        let temp = elem * Gamma_Lanczos((k+1)/2) * HyperGeometricFunction(0.5, 0.5 * (k + 1), 1.5, -Math.pow(elem, 2) / k, k);
        temp /= (Math.sqrt(Math.PI * k) * Gamma_Lanczos(k / 2));
        temp += 0.5;
        //if (temp < y[y.length - 1]) {
        //    alert(y.length);
        //    break;
        //}
        y.push(temp);
    }
    return y;
}
function Student_Density(x, k) {
    let y = [];
    for (let elem of x) {
        let numerator = Gamma_Lanczos((k + 1) / 2);
        let denominator = Math.sqrt(k * Math.PI) * Gamma_Lanczos(k / 2) * Math.pow(1 + elem * elem / k, (k + 1) / 2);
        y.push(numerator / denominator);
    }
    return y;
}
function Fisher_Distribution(d1, d2) {
    let x = [];
    for (let i = 0; i < 20000; i++) {
        do {
            var numerator = Khi_Square_Distribution(d1, 1) / d1;
            var denominator = Khi_Square_Distribution(d2, 1) / d2;
        } while (numerator / denominator > 5)
        x.push(numerator / denominator);
    }
    x.sort((a, b) => a - b);
    return x;
}
function Beta_Function(x, y) {
    return Gamma_Lanczos(x) * Gamma_Lanczos(y) / Gamma_Lanczos(x + y);
}
function Incomplete_Beta_Function(x, a, b) {
    var t1 = Math.pow(x, a);
    var t2 = HyperGeometricFunction(a, 1 - b, a + 1, x, 2);
    var t3 = a;
    return (t1 * t2) / t3;
}
function Regularised_Incomplete_Beta_Function(z, a, b) {
    return Incomplete_Beta_Function(z, a, b) / Beta_Function(a, b);
}
function Fisher_Function(list, d1, d2) {
    let y = [];
    for (let elem of list) {
        y.push(Regularised_Incomplete_Beta_Function(d1 * elem / (d1 * elem + d2), d1 / 2, d2 / 2));
    }
    return y;
}
function Fisher_Density(list, d1, d2) {
    let y = [];
    for (let elem of list) {
        let numerator = Math.sqrt(Math.pow(d1 * elem, d1) * Math.pow(d2, d2) / Math.pow(d1 * elem + d2, d1 + d2));
        let denominator = elem * Beta_Function(d1 / 2, d2 / 2);
        y.push(numerator / denominator);
    }
    return y;
}
function $(id) {
    return document.getElementById(id);
}
function Math_Expectation(distr) {
    if (distr == "uniform") {
        return ($('a').valueAsNumber + $('b').valueAsNumber) / 2;
    }
    else if (distr == "exponential") {
        return Math.pow($('c').valueAsNumber, -1);
    }
    else if (distr == "normal") {
        return $('g').valueAsNumber;
    }
    else if (distr == "khi") {
        return $('d').valueAsNumber;
    }
    else if (distr == "stu") {
        if ($('d').valueAsNumber > 1) {
            return 0;
        }
        return NaN;
    }
    else {
        if ($('e').valueAsNumber > 2) {
            return $('e').valueAsNumber / ($('e').valueAsNumber - 2);
        }
        else return NaN;
    }
}
function Variance(distr) {
    if (distr == "uniform") {
        return Math.pow($('b').valueAsNumber - $('a').valueAsNumber, 2) / 12;
    }
    else if (distr == "exponential") {
        return Math.pow($('c').valueAsNumber, -2);
    }
    else if (distr == "normal") {
        return $('h').valueAsNumber;
    }
    else if (distr == "khi") {
        return 2*$('d').valueAsNumber;
    }
    else if (distr == "stu") {
        if ($('d').valueAsNumber > 2) {
            return $('d').valueAsNumber / ($('d').valueAsNumber-2);
        }
        return NaN;
    }
    else {
        if ($('e').valueAsNumber > 2) {
            return 2 * Math.pow($('e').valueAsNumber, 2) * ($('e').valueAsNumber + $('d').valueAsNumber - 2) / ($('d').valueAsNumber * Math.pow($('e').valueAsNumber - 2, 2) * ($('e').valueAsNumber-4));
        }
        else return NaN;
    }
}

function Deviation(distr) {
    return Math.sqrt(Variance(distr));
}
function Clearcanvas() {
    var canvas = document.getElementById('canvas');
    var c = canvas.getContext('2d');
    c.clearRect(0, 0, canvas.width, canvas.height);
    CreateAxis();
}
function CreateAxis() {
    var canv = document.getElementById("canvas");
    var context = canv.getContext("2d");
    context.beginPath();
    context.strokeStyle = "#000";
    context.font = "bold 14px sans-serif";
    for (let i = 1; i < 10; i++) {
        context.moveTo(0, 500 - i * 50);
        context.fillText(i / 10 + "", 0, 500 - i * 50);
        context.lineTo(1300, 500 - i * 50);
    }
    for (let i = -3; i <= 8; i++) {
        context.moveTo((i + 4) * 100, 500);
        context.lineTo((i + 4) * 100, 0);
        context.fillText(i + "", (i + 4) * 100, 500);
    }
    context.fillText("9", 1290, 500);
    context.fillText("0.0", 0, 500);
    context.fillText("1.0", 0, 13);
    context.stroke();
    context.closePath();
}
function Drawgraphic() {
    if (document.getElementById('1').checked) {
        var x = UniformDistribution(+(document.getElementById('a').value), +(document.getElementById('b').value), 20000);
        if (document.getElementById('f1').checked) {
            var y = UniformDistribution_DistibutionFunction(x, +(document.getElementById('a').value), +(document.getElementById('b').value));

            let c = document.getElementById('canvas').getContext('2d');
            c.beginPath();
            c.lineWidth = "7px";
            c.moveTo(0, 500);
            c.lineTo((x[0] + 4) * 100, 500 - y[0] * 500);
            c.moveTo((x[x.length - 1] + 4) * 100, 0);
            c.lineTo(900, 0);
            c.strokeStyle = "#00f";
            c.stroke();
            c.closePath();
        }
        else var y = Uniform_Density(x, +(document.getElementById('a').value), +(document.getElementById('b').value));
        var me = Math_Expectation('uniform');
        var va = Variance('uniform');
        var dev = Deviation('uniform');
    }
    else if (document.getElementById('2').checked) {
        var x = ExponentialDistribution(+(document.getElementById('c').value), 20000);
        if (document.getElementById('f1').checked) {
            var y = ExponentialDistribution_DistributionFunction(x, +(document.getElementById('c').value));
        }
        else {
            var y = Exp_Density(x, +(document.getElementById('c').value));
        }
        var me = Math_Expectation('exponential');
        var va = Variance('exponential');
        var dev = Deviation('exponential');
    }
    else if (document.getElementById('3').checked) {
        var x = Gauss(20000, $('g').valueAsNumber, $('h').valueAsNumber);
        if (document.getElementById('f1').checked) {
            var y = Gauss_Laplas_Distribution_Distribution_Function(x);
        }
        else var y = Gauss_Density(x);
        var me = Math_Expectation('normal');
        var va = Variance('normal');
        var dev = Deviation('normal');
    }
    else if (document.getElementById('4').checked) {
        var x = Khi_Square_Distribution(document.getElementById('d').valueAsNumber, 20000);
        if (document.getElementById('f1').checked) {
            var y = Khi_Square_Distribution_Distribution_Function(x, document.getElementById('d').valueAsNumber);
        }
        else {
            var y = Khi_Density(x, document.getElementById('d').valueAsNumber);
        }
        var me = Math_Expectation('khi');
        var va = Variance('khi');
        var dev = Deviation('khi');
    }
    else if (document.getElementById('5').checked) {
        var x = Student_Distribution(document.getElementById('d').valueAsNumber);
        if (document.getElementById('f1').checked) {
            var y = Student_Distribution_Distribution_Function(x, document.getElementById('d').valueAsNumber);
        }
        else {
            var y = Student_Density(x, document.getElementById('d').valueAsNumber);
        }
        var me = Math_Expectation('stu');
        var va = Variance('stu');
        var dev = Deviation('stu');
    }
    else if (document.getElementById('6').checked) {
        var x = Fisher_Distribution(document.getElementById('d').valueAsNumber, document.getElementById('e').valueAsNumber);
        if (document.getElementById('f1').checked) {
            var y = Fisher_Function(x, document.getElementById('d').valueAsNumber, document.getElementById('e').valueAsNumber);
        }
        else {
            var y = Fisher_Density(x, document.getElementById('d').valueAsNumber, document.getElementById('e').valueAsNumber);
        }
        var me = Math_Expectation('fish');
        var va = Variance('fish');
        var dev = Deviation('fish');
    }
    else if (document.getElementById('7').checked) {
        var x = CreateZiggurat(20000, $('g').valueAsNumber, $('h').valueAsNumber);
        if (document.getElementById('f1').checked) {
            var y = Gauss_Laplas_Distribution_Distribution_Function(x);
        }
        else {
            var y = Gauss_Density(x);
        }
        var me = Math_Expectation('normal');
        var va = Variance('normal');
        var dev = Deviation('normal');
    }
    for (let el of document.getElementsByClassName('stat')) {
        el.style.display = 'block';
    }
    document.getElementById('me').innerHTML = ("Math expectation: " + me);
    document.getElementById('va').innerHTML = ("Variance: " + va);
    document.getElementById('sig').innerHTML = ("Deviation: " + dev);
    var c = document.getElementById('canvas').getContext('2d');
    c.beginPath();
    c.moveTo((x[0]+4)*100, 500-y[0]*500);
    for (let i = 0; i < x.length; i++) {
        c.lineTo((x[i] + 4) * 100, 500 - y[i] * 500);
    }
    c.strokeStyle = "#00f";
    c.stroke();
    c.closePath();
}
CreateAxis();
