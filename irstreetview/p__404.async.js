(window["webpackJsonp"]=window["webpackJsonp"]||[]).push([[2],{eZYV:function(e,n,t){"use strict";t.d(n,"c",function(){return u}),t.d(n,"a",function(){return r}),t.d(n,"b",function(){return i});var c=t("q1tI");function u(e){var n=Object(c["useRef"])();return Object(c["useEffect"])(function(){n.current=e},[e]),n.current}function r(e,n){void 0===n&&(n=!1);var t=Object(c["useState"])(n),u=t[0],r=t[1];return Object(c["useEffect"])(function(){e?e.style.pointerEvents=u?"none":"all":console.warn("[useFreezeInteraction] element is undefined")},[u]),[u,r]}function i(e,n){var t=Object(c["useRef"])();Object(c["useEffect"])(function(){t.current=e},[e]),Object(c["useEffect"])(function(){function e(){var e;null===(e=t.current)||void 0===e||e.call(t)}if(null!==n){var c=setInterval(e,n);return function(){return clearInterval(c)}}},[n])}},"i6+/":function(e,n,t){"use strict";t.r(n);var c=t("q1tI"),u=t.n(c),r=t("3a4m"),i=t.n(r),f=t("eZYV"),o=function(e){var n=Object(c["useState"])(3),t=n[0],r=n[1];return Object(f["b"])(function(){r(t-1)},1e3),Object(c["useEffect"])(function(){t<=0&&i.a.push("/maps")},[t]),u.a.createElement("div",null,"Unknown route. Redirecting in ",t," seconds",u.a.createElement("div",null,e.children))};n["default"]=o}}]);