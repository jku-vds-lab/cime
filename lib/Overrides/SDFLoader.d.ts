import { Dataset, DatasetType, IVector, Loader, useCancellablePromise } from 'projection-space-explorer';
export declare class SDFLoader implements Loader {
    vectors: IVector[];
    datasetType: DatasetType;
    loadingArea: string;
    resolvePath(entry: any, finished: (dataset: Dataset) => void, cancellablePromise?: ReturnType<typeof useCancellablePromise>['cancellablePromise'], modifiers?: string, controller?: AbortController): void;
    resolveContent(file: any, finished: (dataset: Dataset) => void, cancellablePromise?: ReturnType<typeof useCancellablePromise>['cancellablePromise'], modifiers?: string, controller?: AbortController): void;
    loadCSV(finished: (dataset: Dataset) => void, entry: any, cancellablePromise?: ReturnType<typeof useCancellablePromise>['cancellablePromise'], modifiers?: string, controller?: AbortController, onError?: (error: any) => void): void;
}
//# sourceMappingURL=SDFLoader.d.ts.map