# Decisión de release — mdcheck 1.1.0

## Recomendación técnica

La versión local `1.1.0` es una candidata legítima de release: contiene cambios funcionales delimitados, validación comparativa reproducible, 23 pruebas aprobadas y artefactos que pasan `twine check`.

## Riesgo de publicar sin decisión

Publicarla sin fijar el commit, la fecha del changelog y el vínculo Zenodo produciría una divergencia entre PyPI, CFF, documentación y DOI. Reutilizar el DOI de `1.0.0` sería incorrecto.

## Puerta humana única

Confirmar si `1.1.0` debe convertirse en el siguiente release público. Hasta esa confirmación, el estado correcto es **preparado localmente, no publicado**.
