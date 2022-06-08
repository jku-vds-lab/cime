import { CIMEBackendFromEnv } from "../Backend/CIMEBackend";
import { trackPromise } from "react-promise-tracker";
import {
  AVector,
  CSVLoader,
  Dataset,
  DatasetType,
  IVector,
  Loader,
  useCancellablePromise,
} from "projection-space-explorer";
import * as d3v5 from "d3v5";

function convertFromCSV(vectors) {
  return vectors.map((vector) => {
    return AVector.create(vector);
  });
}

export class SDFLoader implements Loader {
  vectors: IVector[] = [];
  datasetType: DatasetType = DatasetType.None;

  loadingArea = "global_loading_indicator";

  resolvePath(
    entry: any,
    finished: (dataset: Dataset) => void,
    cancellablePromise?: ReturnType<
      typeof useCancellablePromise
    >["cancellablePromise"],
    modifiers?: string,
    controller?: AbortController
  ) {
    if (entry.uploaded) {
      // TODO: Instead of localstorage store it in state?
      // localStorage.setItem("id", entry.path);
      // use file that is already uploaded to backend
      this.loadCSV(finished, entry, cancellablePromise, modifiers, controller);
    } else {
      trackPromise(
        fetch(entry.path, { signal: controller?.signal })
          .then((response) => response.blob())
          .then((result) =>
            this.resolveContent(
              result,
              finished,
              cancellablePromise,
              modifiers,
              controller
            )
          )
          .catch((error) => {
            console.log(error);
          }),
        this.loadingArea
      );
    }
  }

  resolveContent(
    file: any,
    finished: (dataset: Dataset) => void,
    cancellablePromise?: ReturnType<
      typeof useCancellablePromise
    >["cancellablePromise"],
    modifiers?: string,
    controller?: AbortController
  ) {
    const promise = cancellablePromise
      ? cancellablePromise(
          CIMEBackendFromEnv.upload_sdf_file(file, controller),
          controller
        )
      : CIMEBackendFromEnv.upload_sdf_file(file, controller);
    trackPromise(
      promise
        .then((uploaded) => {
          //console.log("UPLOADED", uploaded);
          this.loadCSV(
            finished,
            { display: "", type: this.datasetType, path: uploaded.id },
            cancellablePromise,
            modifiers,
            controller
          );
        })
        .catch((error) => {
          console.log(error);
        }),
      this.loadingArea
    );
  }

  loadCSV(
    finished: (dataset: Dataset) => void,
    entry,
    cancellablePromise?: ReturnType<
      typeof useCancellablePromise
    >["cancellablePromise"],
    modifiers?: string,
    controller?: AbortController,
    onError?: (error: any) => void
  ) {
    // request the server to return a csv file using the unique filename
    const path =
      CIMEBackendFromEnv.baseUrl + "/get_csv/" + entry.path + "/" + modifiers;
    const promise = cancellablePromise
      ? cancellablePromise(
          d3v5.csv(path, {
            ...CIMEBackendFromEnv.fetchParams,
            signal: controller?.signal,
          }),
          controller
        )
      : d3v5.csv(path, {
          ...CIMEBackendFromEnv.fetchParams,
          signal: controller?.signal,
        });
    trackPromise(
      promise
        .then((vectors) => {
          console.log(vectors);
          this.vectors = convertFromCSV(vectors);
          this.datasetType = DatasetType.Chem;
          new CSVLoader()
            .resolve(finished, this.vectors, this.datasetType, entry)
            .catch((error) => {
              onError(error);
            });
        })
        .catch((error) => {
          if (onError) {
            onError(error);
          }
        }),
      this.loadingArea
    );
  }
}
