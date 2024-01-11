---
title: "R-package - To Dos"
format: 
  html:
    self-contained: true
    title-block-banner: true
    page-layout: full
    title-block-banner-color: blue
    keep-md: true

editor: visual
---

::: cell
:::

Last updated: 2024-01-11 by Kenneth

## Package Info

::: panel-tabset
### ToDos

### Documentation

```         
Registered S3 method overwritten by 'printr':
  method                from     
  knit_print.data.frame rmarkdown
```

\<!DOCTYPE html\>

<html>

<head>

<title>R: The 'dsemRLocal' package.</title>

<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />

<meta name="viewport" content="width=device-width, initial-scale=1.0, user-scalable=yes" />

<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.15.3/dist/katex.min.css">

```{=html}
<script type="text/javascript">
const macros = { "\\R": "\\textsf{R}", "\\code": "\\texttt"};
function processMathHTML() {
    var l = document.getElementsByClassName('reqn');
    for (let e of l) { katex.render(e.textContent, e, { throwOnError: false, macros }); }
    return;
}</script>
```
```{=html}
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.15.3/dist/katex.min.js"
    onload="processMathHTML();"></script>
```
<link rel="stylesheet" type="text/css" href="R.css" />

</head>

<body>

::: container
+--------------------+----------------:+
| dsemRLocal-package | R Documentation |
+--------------------+-----------------+

<h2 id="dsemRLocal-package">

The 'dsemRLocal' package.

</h2>

<h3>Description</h3>

<p>A DESCRIPTION OF THE PACKAGE</p>

<h3>References</h3>

<p>Stan Development Team (NA). RStan: the R interface to Stan. R package version 2.32.3. https://mc-stan.org</p>
:::

</body>

</html>
:::

## Workflow

### Modellspezifikation (separiert von den Daten)

Notwendige Schritte:

1.  `VARmodelBuild`

Optionale Schritte:\
Darauf achten, dass alle outputs der Modell-Funktionen miteinander kompatibel sind.

-   `VARmodelConstraints`

-   `VARmodelBetween`

-   `VARmodelMeasurement`

-   `VARmodelPriors`

### Modelschätzung

Eine Funktion zur Modellschätzung. Muss folgendes leisten:

-   Allgemeine Datenvorverarbeitung nach User-Vorgaben

-   Modellspezifische Datenvorverarbeitung (entsprechend den Anforderungen nach Stan-Modelltyp)

-   Aus VARmodel_obj herauslesen, welcher Modelltyp spezifiziert wurde (AR vs. VAR; manifest vs. latent), um entsprechendes Stan-Modell zur Schätzung zu verwenden

## ToDos nach Funktion

### `VARmodelBuild()`

::: {style="float:left"}
<svg width="90.0" height="20" xmlns="http://www.w3.org/2000/svg">

<linearGradient id="a" x2="0" y2="100%"> <stop offset="0" stop-color="#bbb" stop-opacity=".2"/> <stop offset="1" stop-opacity=".1"/> </linearGradient> <rect rx="4" x="0" width="90.0" height="20" fill="#555"/> <rect rx="4" x="0" width="45" height="20" fill="#f0ad4e"/> <rect rx="4" width="90.0" height="20" fill="url(#a)"/> <g fill="#fff" text-anchor="middle" font-family="DejaVu Sans,Verdana,Geneva,sans-serif" font-size="11"> <text x="45.0" y="14"> 50% </text> </g>

</svg>
:::

::: panel-tabset
#### ToDos

<input type="checkbox" checked> Funktion erstellt Datensätze mit Spezifikationen des Structural Models und Measurement Model, falls 1 Indikatoren Funktion packt Datensätze in Liste und speichert diese in einem `VARmodel`-object? </input> <br><input type="checkbox" unchecked> Funktion erstellt Datensätze mit Spezifikationen des Structural Models und Measurement Model, falls 1 Indikatoren Funktion packt Datensätze in Liste und speichert diese in einem `VARmodel`-object?</input> <br>

#### Documentation

\<!DOCTYPE html\>

<html>

<head>

<title>R: Title</title>

<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />

<meta name="viewport" content="width=device-width, initial-scale=1.0, user-scalable=yes" />

<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.15.3/dist/katex.min.css">

```{=html}
<script type="text/javascript">
const macros = { "\\R": "\\textsf{R}", "\\code": "\\texttt"};
function processMathHTML() {
    var l = document.getElementsByClassName('reqn');
    for (let e of l) { katex.render(e.textContent, e, { throwOnError: false, macros }); }
    return;
}</script>
```
```{=html}
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.15.3/dist/katex.min.js"
    onload="processMathHTML();"></script>
```
<link rel="stylesheet" type="text/css" href="R.css" />

</head>

<body>

::: container
+---------------+----------------:+
| VARmodelBuild | R Documentation |
+---------------+-----------------+

<h2 id="VARmodelBuild">

Title

</h2>

<h3>Description</h3>

<p>Title</p>

<h3>Usage</h3>

```{=html}
<pre><code class='language-R'>VARmodelBuild(q, p = NULL)
</code></pre>
```
<h3>Arguments</h3>

|     |                                                                                                                              |
|-----|------------------------------------------------------------------------------------------------------------------------------|
| `q` | integer. The number of time-varying constructs.                                                                              |
| `p` | integer. For multiple-indicator models, specify a vector of length `q` with the number of manifest indicators per construct. |

<h3>Value</h3>

<p>An object of class <code>data.frame</code>.</p>
:::

</body>

</html>
:::

### `VARmodelConstraints()`

::: {style="float:left"}
<svg width="90.0" height="20" xmlns="http://www.w3.org/2000/svg">

<linearGradient id="a" x2="0" y2="100%"> <stop offset="0" stop-color="#bbb" stop-opacity=".2"/> <stop offset="1" stop-opacity=".1"/> </linearGradient> <rect rx="4" x="0" width="90.0" height="20" fill="#555"/> <rect rx="4" x="0" width="45" height="20" fill="#f0ad4e"/> <rect rx="4" width="90.0" height="20" fill="url(#a)"/> <g fill="#fff" text-anchor="middle" font-family="DejaVu Sans,Verdana,Geneva,sans-serif" font-size="11"> <text x="45.0" y="14"> 50% </text> </g>

</svg>
:::

::: panel-tabset
#### ToDos

<input type="checkbox" checked> Funktion erstellt Datensätze mit Spezifikationen des Structural Models und Measurement Model, falls 1 Indikatoren Funktion packt Datensätze in Liste und speichert diese in einem `VARmodel`-object? </input> <br><input type="checkbox" unchecked> Funktion erstellt Datensätze mit Spezifikationen des Structural Models und Measurement Model, falls 1 Indikatoren Funktion packt Datensätze in Liste und speichert diese in einem `VARmodel`-object?</input> <br>

#### Documentation

\<!DOCTYPE html\>

<html>

<head>

<title>R: Title</title>

<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />

<meta name="viewport" content="width=device-width, initial-scale=1.0, user-scalable=yes" />

<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.15.3/dist/katex.min.css">

```{=html}
<script type="text/javascript">
const macros = { "\\R": "\\textsf{R}", "\\code": "\\texttt"};
function processMathHTML() {
    var l = document.getElementsByClassName('reqn');
    for (let e of l) { katex.render(e.textContent, e, { throwOnError: false, macros }); }
    return;
}</script>
```
```{=html}
<script defer src="https://cdn.jsdelivr.net/npm/katex@0.15.3/dist/katex.min.js"
    onload="processMathHTML();"></script>
```
<link rel="stylesheet" type="text/css" href="R.css" />

</head>

<body>

::: container
+---------------------+----------------:+
| VARmodelConstraints | R Documentation |
+---------------------+-----------------+

<h2 id="VARmodelConstraints">

Title

</h2>

<h3>Description</h3>

<p>Title</p>

<h3>Usage</h3>

```{=html}
<pre><code class='language-R'>VARmodelConstraints(
  VARmodel,
  FixInnoVars = F,
  FixInnoCovs = F,
  InnoCovsZero = F,
  FEis0 = NULL,
  REis0 = NULL
)
</code></pre>
```
<h3>Arguments</h3>

|                        |                                                                                                                                                                                                           |
|------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `VARmodel`             | data.frame. Output of VARmodel-Functions.                                                                                                                                                                 |
| `FEis0`                | character. A character vector to indeVARmodel which fiVARmodeled model parameters should be fiVARmodeled to zero (Note: this results in removing the random effect variance of the respective parameter). |
| `REis0`                | logical. Set to TRUE to treat all innovations as independent.                                                                                                                                             |
| `FiVARmodeledInnoVars` | logical. FiVARmodel all random effect variances (eVARmodelcept those of individual traits) to zero.                                                                                                       |
| `FiVARmodeledCovs`     | logical. Set all innovation covariances to a constant value.                                                                                                                                              |
| `CovsZero`             | logical. Set to TRUE to treat all innovations as independent.                                                                                                                                             |

<h3>Value</h3>

<p>An object of class <code>data.frame</code>.</p>
:::

</body>

</html>
:::
