# Estado canónico — mdcheck

Fecha de verificación: 2026-10-12

## Clasificación

Versión local de mantenimiento mayor (`1.1.0`) preparada técnicamente; PyPI todavía sirve `1.0.0`. No está publicada ni debe presentarse como release estable hasta tomar la decisión de versionado.

## Evidencia verificada

- `pyproject.toml` y `CITATION.cff` locales declaran `1.1.0`.
- `CHANGELOG.md` describe los cambios y mantiene la versión como `unreleased`.
- La suite local contiene 23 pruebas aprobadas.
- La orden canónica es `python -m pytest tests -q` ejecutada desde `simcert-suite/mdcheck`; esa ejecución pasó `23 passed`. Ejecutar la misma ruta desde la raíz del contenedor produce errores de importación porque agrupa repositorios independientes.
- Wheel y sdist locales de `1.1.0` pasan `twine check`.
- La comprobación repetida el 2026-10-12 pasó para ambos artefactos: `mdcheck-1.1.0-py3-none-any.whl` y `mdcheck-1.1.0.tar.gz`; ambos contienen metadata de versión 1.1.0.
- Los cambios se apoyan en la validación comparativa de `validation/` frente a pymbar, pyblock y ArviZ.

## Cambios científicos/técnicos de 1.1.0

- Intervalo de confianza corregido para datos autocorrelacionados.
- Consistencia entre réplicas basada en Cochran Q y errores corregidos.
- Parser AMBER compatible con pandas moderno.
- Parser de logs LAMMPS.
- Suite reproducible de validación y búsqueda de equilibración más eficiente.

## Decisión pendiente

Elegir una sola ruta:

1. Publicar `1.1.0`: cambiar `unreleased` por la fecha real, fijar commit/tag, regenerar artefactos y actualizar Zenodo/DOI.
2. Mantener `1.0.0`: conservar la rama local como desarrollo no publicado y no usar sus resultados como evidencia del release público.

No publicar, crear tag ni modificar depósitos externos por inferencia.
